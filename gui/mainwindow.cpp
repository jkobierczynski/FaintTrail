// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "mainwindow.hpp"

#include <QAction>
#include <QApplication>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDir>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QSplitter>
#include <QStatusBar>
#include <QTableWidget>
#include <QThread>
#include <QTimer>
#include <QToolBar>
#include <QVBoxLayout>
#include <QtConcurrent>
#include <cmath>

#include "imageview.hpp"

using namespace fainttrail;

namespace {

enum Layer { LayerFrame = 0, LayerResidual, LayerStatic, LayerSnrMap };

QString sexagesimal(double deg, bool hours) {
    double v = hours ? deg / 15.0 : deg;
    const QChar sign = v < 0 ? QChar('-') : (hours ? QChar() : QChar('+'));
    v = std::fabs(v);
    int a = static_cast<int>(v);
    int b = static_cast<int>((v - a) * 60.0);
    double c = ((v - a) * 60.0 - b) * 60.0;
    if (c >= 59.995) { c = 0; ++b; }
    if (b >= 60) { b = 0; ++a; }
    QString s = QString("%1%2:%3:%4")
                    .arg(sign.isNull() ? QString() : QString(sign))
                    .arg(a, 2, 10, QChar('0'))
                    .arg(b, 2, 10, QChar('0'))
                    .arg(c, hours ? 5 : 4, 'f', hours ? 2 : 1, QChar('0'));
    return s;
}

QTableWidgetItem* numItem(double v, int decimals) {
    auto* it = new QTableWidgetItem;
    it->setData(Qt::DisplayRole, QString::number(v, 'f', decimals).toDouble());
    it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
    return it;
}

QPixmap cutoutPixmap(const Image& cut, int scale) {
    // S/N stacks: fixed range from -2 to 8 sigma so cutouts are comparable
    QImage q = to_qimage(cut, Stretch::Linear, -2.0, 8.0);
    return QPixmap::fromImage(q.scaled(q.width() * scale, q.height() * scale, Qt::IgnoreAspectRatio,
                                       Qt::FastTransformation));
}

}  // namespace

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle(QString("FaintTrail %1").arg(version()));
    buildUi();
    connect(&m_searchWatcher, &QFutureWatcher<SearchOutcome>::finished, this, &MainWindow::searchFinished);
    connect(&m_simWatcher, &QFutureWatcher<SimOutcome>::finished, this, &MainWindow::simulationFinished);
    m_progressTimer = new QTimer(this);
    m_progressTimer->setInterval(100);
    connect(m_progressTimer, &QTimer::timeout, this, &MainWindow::pollProgress);
    m_blinkTimer = new QTimer(this);
    m_blinkTimer->setInterval(220);
    connect(m_blinkTimer, &QTimer::timeout, this, &MainWindow::blinkStep);
    resize(1360, 860);
}

void MainWindow::buildUi() {
    // ---------------------------------------------------------- actions
    auto* fileMenu = menuBar()->addMenu(tr("&File"));
    auto* actOpen = fileMenu->addAction(tr("&Open frames..."), this, &MainWindow::openFiles);
    actOpen->setShortcut(QKeySequence::Open);
    fileMenu->addAction(tr("Open &folder..."), this, &MainWindow::openFolder);
    auto* actSim = fileMenu->addAction(tr("&Simulate sequence..."), this, &MainWindow::simulateDialog);
    m_actImportWcs = fileMenu->addAction(tr("Import &plate solution..."), this, &MainWindow::importPlateSolution);
    m_actImportWcs->setEnabled(false);
    m_actImportWcs->setToolTip(tr("Borrow the WCS of a plate-solved image of the same field (e.g. an ASTAP-solved stack)"));
    fileMenu->addSeparator();
    m_actExportCand = fileMenu->addAction(tr("Export &candidates (CSV)..."), this, &MainWindow::exportCandidates);
    m_actExportTrack = fileMenu->addAction(tr("Export &tracklets (CSV)..."), this, &MainWindow::exportTracklets);
    m_actExportCand->setEnabled(false);
    m_actExportTrack->setEnabled(false);
    fileMenu->addSeparator();
    fileMenu->addAction(tr("&Quit"), qApp, &QApplication::quit)->setShortcut(QKeySequence::Quit);
    menuBar()->addMenu(tr("&Help"))->addAction(tr("&About FaintTrail"), this, &MainWindow::about);

    auto* tb = addToolBar(tr("Main"));
    tb->setMovable(false);
    tb->addAction(actOpen);
    tb->addAction(actSim);
    tb->addAction(m_actImportWcs);

    // ---------------------------------------------------------- parameters dock
    auto* params = new QWidget;
    auto* pl = new QVBoxLayout(params);
    auto* motion = new QGroupBox(tr("Motion range"));
    auto* mf = new QFormLayout(motion);
    m_vmin = new QDoubleSpinBox;
    m_vmin->setRange(0.0, 5000.0);
    m_vmin->setDecimals(1);
    m_vmin->setSuffix(" px/h");
    m_vmin->setValue(SearchParams().vmin_pxh);
    m_vmax = new QDoubleSpinBox;
    m_vmax->setRange(1.0, 5000.0);
    m_vmax->setDecimals(1);
    m_vmax->setSuffix(" px/h");
    m_vmax->setValue(SearchParams().vmax_pxh);
    m_rateLabel = new QLabel;
    m_rateLabel->setWordWrap(true);
    m_rateLabel->setStyleSheet("color: palette(mid);");
    mf->addRow(tr("Slowest"), m_vmin);
    mf->addRow(tr("Fastest"), m_vmax);
    mf->addRow(m_rateLabel);
    connect(m_vmin, qOverload<double>(&QDoubleSpinBox::valueChanged), this, &MainWindow::updateRateLabels);
    connect(m_vmax, qOverload<double>(&QDoubleSpinBox::valueChanged), this, &MainWindow::updateRateLabels);
    pl->addWidget(motion);

    auto* det = new QGroupBox(tr("Detection"));
    auto* df = new QFormLayout(det);
    m_falseRate = new QDoubleSpinBox;
    m_falseRate->setRange(0.0001, 10.0);
    m_falseRate->setDecimals(4);
    m_falseRate->setSingleStep(0.01);
    m_falseRate->setValue(SearchParams().target_false_per_run);
    m_falseRate->setToolTip(tr("Accepted number of pure-noise detections per search. The threshold is "
                               "calibrated on time-shuffled copies of the data to meet it."));
    m_nullRuns = new QSpinBox;
    m_nullRuns->setRange(0, 20);
    m_nullRuns->setValue(SearchParams().null_runs);
    m_nullRuns->setToolTip(tr("Number of time-shuffled null searches used for the calibration."));
    m_threshold = new QDoubleSpinBox;
    m_threshold->setRange(0.0, 100.0);
    m_threshold->setDecimals(2);
    m_threshold->setSpecialValueText(tr("calibrate"));
    m_threshold->setToolTip(tr("Fixed S/N threshold; 'calibrate' derives it from the null searches."));
    m_fwhm = new QDoubleSpinBox;
    m_fwhm->setRange(0.0, 30.0);
    m_fwhm->setDecimals(2);
    m_fwhm->setSuffix(" px");
    m_fwhm->setSpecialValueText(tr("measure"));
    m_threads = new QSpinBox;
    m_threads->setRange(0, 256);
    m_threads->setSpecialValueText(tr("all (%1)").arg(QThread::idealThreadCount()));
    df->addRow(tr("False alarms / search"), m_falseRate);
    df->addRow(tr("Null searches"), m_nullRuns);
    df->addRow(tr("Fixed threshold"), m_threshold);
    df->addRow(tr("PSF FWHM"), m_fwhm);
    df->addRow(tr("Threads"), m_threads);
    m_backend = new QComboBox;
    m_backend->addItem(tr("auto"), QString("auto"));
    for (const auto& n : backend_names()) {
        const QString q = QString::fromStdString(n);
        m_backend->addItem(q, q);
        m_backend->setItemData(m_backend->count() - 1, QString::fromStdString(make_backend(n)->description()), Qt::ToolTipRole);
    }
    m_backend->setToolTip(tr("Implementation of the stacking loop"));
    df->addRow(tr("Backend"), m_backend);
    pl->addWidget(det);

    auto* runRow = new QHBoxLayout;
    m_searchButton = new QPushButton(tr("Search"));
    m_searchButton->setDefault(true);
    m_searchButton->setEnabled(false);
    m_cancelButton = new QPushButton(tr("Cancel"));
    m_cancelButton->setEnabled(false);
    connect(m_searchButton, &QPushButton::clicked, this, &MainWindow::runSearch);
    connect(m_cancelButton, &QPushButton::clicked, this, &MainWindow::cancelSearch);
    runRow->addWidget(m_searchButton);
    runRow->addWidget(m_cancelButton);
    pl->addLayout(runRow);

    m_summary = new QLabel(tr("No data loaded."));
    m_summary->setWordWrap(true);
    m_summary->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_summary->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    pl->addWidget(m_summary);
    pl->addStretch(1);

    auto* paramDock = new QDockWidget(tr("Search"), this);
    paramDock->setObjectName("search");
    paramDock->setFeatures(QDockWidget::DockWidgetMovable);
    paramDock->setWidget(params);
    params->setMinimumWidth(250);
    addDockWidget(Qt::LeftDockWidgetArea, paramDock);

    // ---------------------------------------------------------- candidate detail dock
    auto* detail = new QWidget;
    auto* dl = new QVBoxLayout(detail);
    auto* cutRow = new QHBoxLayout;
    const char* cutTitles[3] = {"all frames", "first half", "second half"};
    for (int i = 0; i < 3; ++i) {
        auto* col = new QVBoxLayout;
        m_cut[i] = new QLabel;
        m_cut[i]->setFixedSize(105, 105);
        m_cut[i]->setAlignment(Qt::AlignCenter);
        m_cut[i]->setStyleSheet("background: #121418; border: 1px solid #333;");
        auto* t = new QLabel(tr(cutTitles[i]));
        t->setAlignment(Qt::AlignCenter);
        col->addWidget(m_cut[i]);
        col->addWidget(t);
        cutRow->addLayout(col);
    }
    dl->addLayout(cutRow);
    m_details = new QLabel(tr("Select a candidate to inspect its stacks."));
    m_details->setWordWrap(true);
    m_details->setTextFormat(Qt::RichText);
    m_details->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_details->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    dl->addWidget(m_details);
    dl->addStretch(1);
    auto* detailDock = new QDockWidget(tr("Candidate"), this);
    detailDock->setObjectName("candidate");
    detailDock->setFeatures(QDockWidget::DockWidgetMovable);
    detailDock->setWidget(detail);
    addDockWidget(Qt::RightDockWidgetArea, detailDock);

    // ---------------------------------------------------------- central: image + table
    auto* center = new QWidget;
    auto* cl = new QVBoxLayout(center);
    cl->setContentsMargins(4, 4, 4, 4);
    auto* viewBar = new QHBoxLayout;
    m_layer = new QComboBox;
    m_layer->addItem(tr("Frame"));
    m_layer->addItem(tr("Residual (moving and changing only)"));
    m_layer->addItem(tr("Static sky"));
    m_layer->addItem(tr("Max S/N over all motions"));
    connect(m_layer, qOverload<int>(&QComboBox::currentIndexChanged), this, &MainWindow::layerChanged);
    m_frameSlider = new QSlider(Qt::Horizontal);
    m_frameSlider->setEnabled(false);
    connect(m_frameSlider, &QSlider::valueChanged, this, &MainWindow::frameChanged);
    m_frameLabel = new QLabel("-");
    m_frameLabel->setMinimumWidth(190);
    m_blinkButton = new QPushButton(tr("Blink"));
    m_blinkButton->setCheckable(true);
    m_blinkButton->setEnabled(false);
    connect(m_blinkButton, &QPushButton::toggled, this, &MainWindow::toggleBlink);
    viewBar->addWidget(new QLabel(tr("Show")));
    viewBar->addWidget(m_layer);
    viewBar->addSpacing(12);
    viewBar->addWidget(m_frameSlider, 1);
    viewBar->addWidget(m_frameLabel);
    viewBar->addWidget(m_blinkButton);
    cl->addLayout(viewBar);

    auto* split = new QSplitter(Qt::Vertical);
    m_view = new ImageView;
    connect(m_view, &ImageView::hovered, this, [this](double x, double y, double v, bool inside) {
        if (!inside) {
            m_hover->clear();
            return;
        }
        QString s = QString("x %1  y %2  value %3").arg(x + 1.0, 0, 'f', 1).arg(y + 1.0, 0, 'f', 1).arg(v, 0, 'g', 5);
        Wcs w;
        if (m_prep) w = m_prep->wcs;
        else if (!m_frames.empty()) w = m_frames.front().wcs;
        double ra, dec;
        if (w.pix2sky(x, y, ra, dec)) s += QString("   RA %1  Dec %2").arg(sexagesimal(ra, true), sexagesimal(dec, false));
        m_hover->setText(s);
    });
    connect(m_view, &ImageView::markerClicked, this, [this](int i) {
        // markers are candidates first (in table order), truth markers after
        if (m_result && i < static_cast<int>(m_result->candidates.size())) selectCandidate(i);
    });
    split->addWidget(m_view);
    m_table = new QTableWidget;
    const QStringList headers = {tr("#"),  tr("pass"), tr("x"),      tr("y"),     tr("rate \"/min"), tr("PA °"),
                                 tr("S/N"), tr("half 1"), tr("half 2"), tr("max frame"), tr("FA/search"), tr("mag"),
                                 tr("flags")};
    m_table->setColumnCount(headers.size());
    m_table->setHorizontalHeaderLabels(headers);
    m_table->verticalHeader()->setVisible(false);
    m_table->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_table->setSelectionMode(QAbstractItemView::SingleSelection);
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->horizontalHeader()->setStretchLastSection(true);
    m_table->setAlternatingRowColors(true);
    connect(m_table, &QTableWidget::itemSelectionChanged, this, &MainWindow::tableSelectionChanged);
    split->addWidget(m_table);
    m_table->setMinimumHeight(150);
    split->setStretchFactor(0, 3);
    split->setStretchFactor(1, 1);
    cl->addWidget(split, 1);
    setCentralWidget(center);

    // ---------------------------------------------------------- status bar
    m_status = new QLabel(tr("Ready"));
    m_hover = new QLabel;
    m_progressBar = new QProgressBar;
    m_progressBar->setRange(0, 1000);
    m_progressBar->setMaximumWidth(240);
    m_progressBar->setVisible(false);
    statusBar()->addWidget(m_status, 1);
    statusBar()->addPermanentWidget(m_hover);
    statusBar()->addPermanentWidget(m_progressBar);
    updateRateLabels();
}

// ------------------------------------------------------------------ data in
void MainWindow::openFiles() {
    const QStringList files = QFileDialog::getOpenFileNames(this, tr("Open registered frames"), QString(),
                                                            tr("FITS images (*.fits *.fit *.fts);;All files (*)"));
    if (!files.isEmpty()) loadFiles(files);
}

void MainWindow::openFolder() {
    const QString dir = QFileDialog::getExistingDirectory(this, tr("Folder with registered FITS frames"));
    if (dir.isEmpty()) return;
    QStringList files;
    for (const auto& fi : QDir(dir).entryInfoList({"*.fits", "*.fit", "*.fts", "*.FITS", "*.FIT", "*.FTS"}, QDir::Files,
                                                  QDir::Name))
        files << fi.absoluteFilePath();
    if (files.isEmpty()) {
        QMessageBox::information(this, tr("No frames"), tr("No FITS files found in %1.").arg(dir));
        return;
    }
    loadFiles(files);
}

bool MainWindow::loadFiles(const QStringList& files) {
    std::vector<std::string> paths;
    for (const auto& f : files) paths.push_back(f.toStdString());
    QApplication::setOverrideCursor(Qt::WaitCursor);
    try {
        auto frames = load_frames(paths);
        QApplication::restoreOverrideCursor();
        m_truth.clear();
        setFrames(std::move(frames), QFileInfo(files.front()).absolutePath());
        return true;
    } catch (const std::exception& e) {
        QApplication::restoreOverrideCursor();
        QMessageBox::warning(this, tr("Cannot load frames"), QString::fromStdString(e.what()));
        return false;
    }
}

void MainWindow::simulateDialog() {
    QDialog dlg(this);
    dlg.setWindowTitle(tr("Simulate an observing sequence"));
    auto* form = new QFormLayout(&dlg);
    SimParams d;
    auto spin = [&](int lo, int hi, int v) {
        auto* s = new QSpinBox;
        s->setRange(lo, hi);
        s->setValue(v);
        return s;
    };
    auto dspin = [&](double lo, double hi, double v, int dec, const QString& suffix) {
        auto* s = new QDoubleSpinBox;
        s->setRange(lo, hi);
        s->setDecimals(dec);
        s->setValue(v);
        s->setSuffix(suffix);
        return s;
    };
    auto* w = spin(64, 4096, d.width);
    auto* h = spin(64, 4096, d.height);
    auto* n = spin(5, 500, d.nframes);
    auto* exp = dspin(1, 600, d.exptime_s, 1, " s");
    auto* cad = dspin(1, 600, d.cadence_s, 1, " s");
    auto* fwhm = dspin(1, 15, d.fwhm_px, 2, " px");
    auto* scale = dspin(0.1, 20, d.pixscale_arcsec, 2, " \"/px");
    auto* stars = spin(0, 100000, d.nstars);
    auto* movers = spin(0, 200, d.nmovers);
    auto* magB = dspin(10, 25, d.mover_mag_bright, 2, "");
    auto* magF = dspin(10, 25, d.mover_mag_faint, 2, "");
    auto* rmin = dspin(0.01, 100, d.rate_min_arcsec_min, 2, " \"/min");
    auto* rmax = dspin(0.01, 100, d.rate_max_arcsec_min, 2, " \"/min");
    auto* seed = spin(1, 1000000, static_cast<int>(d.seed));
    form->addRow(tr("Width"), w);
    form->addRow(tr("Height"), h);
    form->addRow(tr("Frames"), n);
    form->addRow(tr("Exposure"), exp);
    form->addRow(tr("Cadence"), cad);
    form->addRow(tr("PSF FWHM"), fwhm);
    form->addRow(tr("Pixel scale"), scale);
    form->addRow(tr("Stars"), stars);
    form->addRow(tr("Moving objects"), movers);
    form->addRow(tr("Brightest mover (mag)"), magB);
    form->addRow(tr("Faintest mover (mag)"), magF);
    form->addRow(tr("Slowest rate"), rmin);
    form->addRow(tr("Fastest rate"), rmax);
    form->addRow(tr("Random seed"), seed);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    buttons->button(QDialogButtonBox::Ok)->setText(tr("Simulate and search"));
    connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);
    form->addRow(buttons);
    if (dlg.exec() != QDialog::Accepted) return;
    if (cad->value() < exp->value()) cad->setValue(exp->value());
    SimParams p;
    p.width = w->value();
    p.height = h->value();
    p.nframes = n->value();
    p.exptime_s = exp->value();
    p.cadence_s = cad->value();
    p.fwhm_px = fwhm->value();
    p.pixscale_arcsec = scale->value();
    p.nstars = stars->value();
    p.nmovers = movers->value();
    p.mover_mag_bright = std::min(magB->value(), magF->value());
    p.mover_mag_faint = std::max(magB->value(), magF->value());
    p.rate_min_arcsec_min = std::min(rmin->value(), rmax->value());
    p.rate_max_arcsec_min = std::max(rmin->value(), rmax->value());
    p.seed = static_cast<std::uint64_t>(seed->value());
    simulate(p, true);
}

void MainWindow::simulate(const SimParams& params, bool then_search) {
    if (m_searchWatcher.isRunning() || m_simWatcher.isRunning()) return;
    m_searchAfterSim = then_search;
    setBusy(true, tr("Simulating %1 frames...").arg(params.nframes));
    m_simWatcher.setFuture(QtConcurrent::run([params]() {
        SimOutcome out;
        out.params = params;
        try {
            out.sim = std::make_shared<SimResult>(fainttrail::simulate(params));
        } catch (const std::exception& e) {
            out.error = QString::fromStdString(e.what());
        }
        return out;
    }));
}

void MainWindow::simulationFinished() {
    setBusy(false, QString());
    const SimOutcome out = m_simWatcher.result();
    if (!out.sim) {
        QMessageBox::warning(this, tr("Simulation failed"), out.error);
        return;
    }
    m_truth = out.sim->truth;
    // match the search range to the simulated rates
    const double scale = out.params.pixscale_arcsec;
    m_vmin->setValue(std::max(0.0, std::floor(arcsec_per_min_to_px_per_h(out.params.rate_min_arcsec_min, scale) * 0.9)));
    m_vmax->setValue(std::ceil(arcsec_per_min_to_px_per_h(out.params.rate_max_arcsec_min, scale) * 1.1));
    setFrames(std::move(out.sim->frames), tr("simulation (seed %1)").arg(out.params.seed));
    if (m_searchAfterSim) runSearch();
}

void MainWindow::setFrames(std::vector<Frame> frames, const QString& origin) {
    m_blinkButton->setChecked(false);
    m_frames = std::move(frames);
    clearResults();
    const int n = frameCount();
    m_frameSlider->blockSignals(true);
    m_frameSlider->setRange(0, std::max(0, n - 1));
    m_frameSlider->setValue(n / 2);
    m_frameSlider->blockSignals(false);
    m_frameSlider->setEnabled(n > 0);
    m_blinkButton->setEnabled(n > 1);
    m_searchButton->setEnabled(n >= 5);
    m_actImportWcs->setEnabled(n > 0);
    m_layer->setCurrentIndex(LayerFrame);
    if (n > 0) {
        const double span_min = (m_frames.back().mjd_mid - m_frames.front().mjd_mid) * 1440.0;
        QString s = tr("<b>%1 frames</b> of %2 × %3 px from %4<br>span %5 min, %6 UTC to %7 UTC")
                        .arg(n)
                        .arg(m_frames[0].image.width)
                        .arg(m_frames[0].image.height)
                        .arg(origin.toHtmlEscaped())
                        .arg(span_min, 0, 'f', 1)
                        .arg(QString::fromStdString(iso_from_mjd(m_frames.front().mjd_mid, 0)).mid(11))
                        .arg(QString::fromStdString(iso_from_mjd(m_frames.back().mjd_mid, 0)).mid(11));
        if (pixelScale() > 0) s += tr("<br>pixel scale %1\"/px (WCS)").arg(pixelScale(), 0, 'f', 2);
        if (!m_truth.empty()) s += tr("<br>%1 simulated movers shown as blue crosses").arg(m_truth.size());
        m_summary->setText(s);
    }
    m_view->fitToWindow();
    refreshImage();
    refreshMarkers();
    updateRateLabels();
    m_status->setText(tr("Loaded %1 frames").arg(n));
}

void MainWindow::clearResults() {
    m_prep.reset();
    m_result.reset();
    m_selected = -1;
    m_table->setRowCount(0);
    for (auto* c : m_cut) c->clear();
    m_details->setText(tr("Select a candidate to inspect its stacks."));
    m_actExportCand->setEnabled(false);
    m_actExportTrack->setEnabled(false);
}

// ------------------------------------------------------------------ search
SearchParams MainWindow::currentParams() const {
    SearchParams p;
    p.vmin_pxh = std::min(m_vmin->value(), m_vmax->value());
    p.vmax_pxh = std::max(m_vmin->value(), m_vmax->value());
    p.target_false_per_run = m_falseRate->value();
    p.null_runs = m_nullRuns->value();
    p.threshold_snr = m_threshold->value();
    p.fwhm_px = m_fwhm->value();
    p.threads = m_threads->value();
    p.backend = m_backend->currentData().toString().toStdString();
    return p;
}

void MainWindow::runSearch() {
    if (frameCount() < 5 || m_searchWatcher.isRunning()) return;
    m_blinkButton->setChecked(false);
    const SearchParams params = currentParams();
    if (params.threshold_snr <= 0.0 && params.null_runs == 0 &&
        QMessageBox::question(this, tr("No calibration"),
                              tr("Without null searches the threshold comes from a Gaussian approximation, which "
                                 "underestimates real false alarms. Continue?")) != QMessageBox::Yes)
        return;
    clearResults();
    m_cancel = false;
    {
        QMutexLocker lock(&m_progressMutex);
        m_progress = 0.0;
        m_stage = tr("starting");
    }
    setBusy(true, tr("Searching..."));
    auto frames = std::make_shared<std::vector<Frame>>(m_frames);
    m_searchWatcher.setFuture(QtConcurrent::run([this, frames, params]() {
        SearchOutcome out;
        auto cb = [this](double f, const std::string& stage) {
            {
                QMutexLocker lock(&m_progressMutex);
                m_progress = f;
                m_stage = QString::fromStdString(stage);
            }
            return !m_cancel.load();
        };
        try {
            out.prep = std::make_shared<Prepared>(
                prepare(*frames, params, [&](double f, const std::string& s) { return cb(0.1 * f, s); }));
            out.result = std::make_shared<SearchResult>(
                search(*out.prep, params, [&](double f, const std::string& s) { return cb(0.1 + 0.9 * f, s); }));
        } catch (const Cancelled&) {
            out.cancelled = true;
        } catch (const std::exception& e) {
            out.error = QString::fromStdString(e.what());
        }
        return out;
    }));
}

void MainWindow::cancelSearch() {
    m_cancel = true;
    m_status->setText(tr("Cancelling..."));
}

void MainWindow::pollProgress() {
    QMutexLocker lock(&m_progressMutex);
    m_progressBar->setValue(static_cast<int>(m_progress * 1000));
    if (m_searchWatcher.isRunning()) m_status->setText(tr("Searching: %1").arg(m_stage));
}

void MainWindow::searchFinished() {
    setBusy(false, QString());
    const SearchOutcome out = m_searchWatcher.result();
    if (out.cancelled) {
        m_status->setText(tr("Search cancelled"));
        return;
    }
    if (!out.result) {
        m_status->setText(tr("Search failed"));
        QMessageBox::warning(this, tr("Search failed"), out.error);
        return;
    }
    m_prep = out.prep;
    m_result = out.result;
    const auto& r = *m_result;
    int ok = 0;
    for (const auto& c : r.candidates) ok += c.ok;
    m_status->setText(tr("%1 candidates passed, %2 flagged, %3 s").arg(ok).arg(r.candidates.size() - ok).arg(r.elapsed_s, 0, 'f', 1));

    QString s = m_summary->text();
    const int cut = s.indexOf("<hr>");
    if (cut >= 0) s = s.left(cut);
    s += tr("<hr><b>Search</b><br>%1 trial motions, grid step %2 px/h<br>PSF FWHM %3 px<br>"
            "threshold S/N <b>%4</b> (%5)<br>slowest reliable motion %6 px/h<br>backend %9<br>"
            "<b>%7 passed</b>, %8 flagged")
             .arg(r.n_velocities)
             .arg(r.dv_pxh, 0, 'f', 2)
             .arg(r.fwhm_px, 0, 'f', 2)
             .arg(r.threshold_snr, 0, 'f', 2)
             .arg(QString::fromStdString(r.threshold_method))
             .arg(r.slowest_reliable_pxh, 0, 'f', 1)
             .arg(ok)
             .arg(r.candidates.size() - ok)
             .arg(QString::fromStdString(r.backend).toHtmlEscaped());
    if (!m_truth.empty()) {
        std::vector<SimMover> movers;
        for (const auto& t : m_truth) movers.push_back(t.mover);
        const auto m = match(r.candidates, movers, 2.0, std::max(3.0, 2.0 * r.dv_pxh), true);
        s += tr("<br>truth: %1 of %2 movers recovered, %3 passing candidates unmatched")
                 .arg(m.size())
                 .arg(m_truth.size())
                 .arg(ok - static_cast<int>(m.size()));
    }
    m_summary->setText(s);
    m_actExportCand->setEnabled(true);
    m_actExportTrack->setEnabled(true);
    refreshTable();
    refreshImage();
    refreshMarkers();
    if (!r.candidates.empty()) selectCandidate(0);
    if (m_onSearchFinished) m_onSearchFinished();
}

void MainWindow::setBusy(bool busy, const QString& what) {
    m_searchButton->setEnabled(!busy && frameCount() >= 5);
    m_cancelButton->setEnabled(busy && m_searchWatcher.isRunning());
    m_progressBar->setVisible(busy);
    m_progressBar->setValue(0);
    if (busy) {
        m_status->setText(what);
        m_progressTimer->start();
    } else {
        m_progressTimer->stop();
    }
}

// ------------------------------------------------------------------ display
double MainWindow::pixelScale() const {
    if (m_prep && m_prep->wcs.valid) return m_prep->wcs.pixel_scale_arcsec();
    for (const auto& f : m_frames)
        if (f.wcs.valid) return f.wcs.pixel_scale_arcsec();
    return 0.0;
}

void MainWindow::updateRateLabels() {
    const double s = pixelScale();
    if (s > 0)
        m_rateLabel->setText(tr("= %1 to %2 \"/min at %3\"/px")
                                 .arg(px_per_h_to_arcsec_per_min(m_vmin->value(), s), 0, 'f', 3)
                                 .arg(px_per_h_to_arcsec_per_min(m_vmax->value(), s), 0, 'f', 3)
                                 .arg(s, 0, 'f', 2));
    else
        m_rateLabel->setText(tr("No WCS: rates in pixels per hour"));
}

void MainWindow::layerChanged() {
    const int layer = m_layer->currentIndex();
    if ((layer == LayerResidual || layer == LayerStatic) && !m_prep) {
        m_status->setText(tr("Run a search first to see this layer"));
    } else if (layer == LayerSnrMap && !m_result) {
        m_status->setText(tr("Run a search first to see this layer"));
    }
    const bool per_frame = layer == LayerFrame || layer == LayerResidual;
    m_frameSlider->setEnabled(per_frame && frameCount() > 0);
    m_blinkButton->setEnabled(per_frame && frameCount() > 1);
    refreshImage();
    refreshMarkers();
}

void MainWindow::frameChanged(int) {
    refreshImage();
    refreshMarkers();
}

void MainWindow::toggleBlink() {
    if (m_blinkButton->isChecked()) m_blinkTimer->start();
    else m_blinkTimer->stop();
}

void MainWindow::blinkStep() {
    if (frameCount() < 2) return;
    m_frameSlider->setValue((m_frameSlider->value() + 1) % frameCount());
}

void MainWindow::refreshImage() {
    if (m_frames.empty()) return;
    const int layer = m_layer->currentIndex();
    const int i = std::clamp(m_frameSlider->value(), 0, frameCount() - 1);
    const Frame& f = m_frames[static_cast<std::size_t>(i)];
    m_frameLabel->setText(tr("frame %1/%2  %3").arg(i + 1).arg(frameCount()).arg(
        QString::fromStdString(iso_from_mjd(f.mjd_mid, 1)).mid(11)));
    switch (layer) {
        case LayerResidual:
            if (m_prep) {
                const double sg = m_prep->sigma[static_cast<std::size_t>(i)];
                m_view->setImage(m_prep->residual[static_cast<std::size_t>(i)], Stretch::Linear, -3 * sg, 6 * sg);
                return;
            }
            break;
        case LayerStatic:
            if (m_prep) {
                m_view->setImage(m_prep->static_sky, Stretch::Asinh);
                return;
            }
            break;
        case LayerSnrMap:
            if (m_result) {
                m_view->setImage(m_result->max_snr, Stretch::Linear, 2.0, std::max(8.0, m_result->threshold_snr + 2));
                return;
            }
            break;
        default:
            break;
    }
    m_view->setImage(f.image, Stretch::Asinh);
}

void MainWindow::refreshMarkers() {
    std::vector<Marker> markers;
    const int layer = m_layer->currentIndex();
    const bool per_frame = layer == LayerFrame || layer == LayerResidual;
    double t = 0.0;  // hours from the reference epoch
    double t_ref = 0.0;
    if (m_result) t_ref = m_result->t_ref_mjd;
    else if (!m_frames.empty()) t_ref = 0.5 * (m_frames.front().mjd_mid + m_frames.back().mjd_mid);
    if (per_frame && !m_frames.empty())
        t = (m_frames[static_cast<std::size_t>(std::clamp(m_frameSlider->value(), 0, frameCount() - 1))].mjd_mid - t_ref) * 24.0;
    if (m_result) {
        const double r = m_result->fwhm_px * 2.0;
        for (std::size_t k = 0; k < m_result->candidates.size(); ++k) {
            const auto& c = m_result->candidates[k];
            Marker m;
            m.x = c.x + c.vx_pxh * t;
            m.y = c.y + c.vy_pxh * t;
            m.radius = r;
            m.color = c.ok ? QColor(80, 220, 120) : QColor(240, 160, 60);
            m.label = QString::number(c.id);
            m.selected = static_cast<int>(k) == m_selected;
            if (m.selected && m_prep) {
                const double t0 = m_prep->t_h.front(), t1 = m_prep->t_h.back();
                m.track = {QPointF(c.x + c.vx_pxh * t0, c.y + c.vy_pxh * t0), QPointF(c.x + c.vx_pxh * t1, c.y + c.vy_pxh * t1)};
            }
            if (!c.ok && !m.selected) m.color.setAlpha(170);
            markers.push_back(m);
        }
    }
    for (const auto& tr_ : m_truth) {
        Marker m;
        m.x = tr_.mover.x_ref + tr_.mover.vx_pxh * t;
        m.y = tr_.mover.y_ref + tr_.mover.vy_pxh * t;
        m.radius = 4;
        m.cross = true;
        m.color = QColor(90, 160, 255);
        markers.push_back(m);
    }
    m_view->setMarkers(std::move(markers));
}

void MainWindow::refreshTable() {
    m_table->setSortingEnabled(false);
    m_table->setRowCount(0);
    if (!m_result) return;
    const auto& cs = m_result->candidates;
    m_table->setRowCount(static_cast<int>(cs.size()));
    const bool rates = pixelScale() > 0;
    for (int row = 0; row < static_cast<int>(cs.size()); ++row) {
        const auto& c = cs[static_cast<std::size_t>(row)];
        int col = 0;
        m_table->setItem(row, col++, numItem(c.id, 0));
        auto* okItem = new QTableWidgetItem(c.ok ? QString::fromUtf8("✓") : QString::fromUtf8("✗"));
        okItem->setTextAlignment(Qt::AlignCenter);
        okItem->setForeground(c.ok ? QColor(60, 180, 90) : QColor(220, 130, 40));
        m_table->setItem(row, col++, okItem);
        m_table->setItem(row, col++, numItem(c.x + 1.0, 2));
        m_table->setItem(row, col++, numItem(c.y + 1.0, 2));
        m_table->setItem(row, col++, numItem(rates ? c.rate_arcsec_min : std::hypot(c.vx_pxh, c.vy_pxh), 3));
        if (rates) m_table->setItem(row, col++, numItem(c.pa_deg, 1));
        else m_table->setItem(row, col++, new QTableWidgetItem(QString::fromUtf8("\u2013")));
        m_table->setItem(row, col++, numItem(c.snr, 2));
        m_table->setItem(row, col++, numItem(c.snr_half1, 2));
        m_table->setItem(row, col++, numItem(c.snr_half2, 2));
        m_table->setItem(row, col++, numItem(c.max_frame_fraction, 3));
        m_table->setItem(row, col++, new QTableWidgetItem(QString::number(c.false_alarms_per_run, 'g', 2)));
        m_table->setItem(row, col++, numItem(c.mag_inst, 2));
        m_table->setItem(row, col++, new QTableWidgetItem(QString::fromStdString(c.flags)));
    }
    if (!rates) m_table->horizontalHeaderItem(4)->setText(tr("v px/h"));
    m_table->resizeColumnsToContents();
}

void MainWindow::tableSelectionChanged() {
    const auto rows = m_table->selectionModel()->selectedRows();
    if (rows.isEmpty()) return;
    const int row = rows.front().row();
    if (row == m_selected) return;
    m_selected = row;
    refreshMarkers();
    refreshDetails();
    if (m_result && row >= 0 && row < static_cast<int>(m_result->candidates.size())) {
        const auto& c = m_result->candidates[static_cast<std::size_t>(row)];
        const int layer = m_layer->currentIndex();
        double t = 0.0;
        if ((layer == LayerFrame || layer == LayerResidual) && m_prep)
            t = m_prep->t_h[static_cast<std::size_t>(std::clamp(m_frameSlider->value(), 0, frameCount() - 1))];
        m_view->centerOn(c.x + c.vx_pxh * t, c.y + c.vy_pxh * t);
    }
}

void MainWindow::selectCandidate(int row) {
    if (!m_result || row < 0 || row >= m_table->rowCount()) return;
    m_table->selectRow(row);
    m_table->scrollToItem(m_table->item(row, 0));
}

void MainWindow::refreshDetails() {
    if (!m_result || !m_prep || m_selected < 0 || m_selected >= static_cast<int>(m_result->candidates.size())) return;
    const auto& c = m_result->candidates[static_cast<std::size_t>(m_selected)];
    for (int k = 0; k < 3; ++k) {
        const Image cut = stack_cutout(*m_prep, c.x, c.y, c.vx_pxh, c.vy_pxh, 10, true, k);
        m_cut[k]->setPixmap(cutoutPixmap(cut, 5));
    }
    const double scale = pixelScale();
    QString html = tr("<b>Candidate %1</b> — %2<br>").arg(c.id).arg(
        c.ok ? tr("<span style='color:#3cb45a'>passed all tests</span>")
             : tr("<span style='color:#dc8228'>flagged: %1</span>").arg(QString::fromStdString(c.flags)));
    html += tr("reference epoch %1 UTC<br>").arg(QString::fromStdString(iso_from_mjd(m_result->t_ref_mjd, 2)));
    html += tr("position x %1, y %2 (± %3 px)<br>").arg(c.x + 1.0, 0, 'f', 2).arg(c.y + 1.0, 0, 'f', 2).arg(c.pos_err_px, 0, 'f', 2);
    if (m_prep->wcs.valid)
        html += tr("RA %1 Dec %2 (± %3\")<br>").arg(sexagesimal(c.ra_deg, true), sexagesimal(c.dec_deg, false)).arg(c.pos_err_px * scale, 0, 'f', 2);
    html += tr("motion %1, %2 px/h (± %3)<br>").arg(c.vx_pxh, 0, 'f', 2).arg(c.vy_pxh, 0, 'f', 2).arg(c.vel_err_pxh, 0, 'f', 2);
    if (scale > 0) html += tr("rate %1\"/min, PA %2°<br>").arg(c.rate_arcsec_min, 0, 'f', 3).arg(c.pa_deg, 0, 'f', 1);
    html += tr("<br>S/N %1 (grid %2, threshold %3)<br>").arg(c.snr, 0, 'f', 2).arg(c.snr_coarse, 0, 'f', 2).arg(m_result->threshold_snr, 0, 'f', 2);
    html += tr("halves %1 / %2 (expected ≈ %3 each)<br>").arg(c.snr_half1, 0, 'f', 2).arg(c.snr_half2, 0, 'f', 2).arg(c.snr_coarse / std::sqrt(2.0), 0, 'f', 2);
    html += tr("largest single-frame share %1<br>").arg(c.max_frame_fraction, 0, 'f', 3);
    html += tr("expected noise peaks this strong: %1 per search<br>").arg(c.false_alarms_per_run, 0, 'g', 2);
    html += tr("instrumental mag %1 (flux %2)").arg(c.mag_inst, 0, 'f', 2).arg(c.flux, 0, 'f', 1);
    if (!c.flags.empty()) {
        html += tr("<br><br><small>half: one half of the sequence lacks the signal · single: one frame dominates · "
                   "ghost: sidelobe of a brighter object · edge: track runs along the blank border of registered frames · detector: fixed on the sensor (warm pixel) · masked: track mostly over masked pixels · "
                   "vedge: at the fastest searched motion · slow: below the reliable motion limit</small>");
    }
    m_details->setText(html);
}

// ------------------------------------------------------------------ plate solution
void MainWindow::importPlateSolution() {
    if (m_frames.empty()) return;
    const QString path = QFileDialog::getOpenFileName(this, tr("Plate-solved image of the same field"), QString(),
                                                      tr("FITS images (*.fits *.fit *.fts *.wcs);;All files (*)"));
    if (!path.isEmpty()) importPlateSolutionFrom(path);
}

bool MainWindow::importPlateSolutionFrom(const QString& path) {
    if (m_frames.empty()) return false;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    try {
        const TransferReport rep = apply_wcs_from_file(path.toStdString(), m_frames);
        QApplication::restoreOverrideCursor();
        clearResults();
        QString s = m_summary->text();
        const int cut = s.indexOf("<hr>");
        if (cut >= 0) s = s.left(cut);
        s += tr("<br>plate solution from %1: %2 stars, rms %3 px, %4\"/px")
                 .arg(QFileInfo(path).fileName().toHtmlEscaped())
                 .arg(rep.matched)
                 .arg(rep.rms_px, 0, 'f', 2)
                 .arg(m_frames.front().wcs.pixel_scale_arcsec(), 0, 'f', 3);
        if (rep.source_has_distortion) s += tr(" (distortion terms ignored)");
        m_summary->setText(s);
        updateRateLabels();
        m_status->setText(tr("Plate solution imported; search again for RA/Dec and rates"));
        return true;
    } catch (const std::exception& e) {
        QApplication::restoreOverrideCursor();
        QMessageBox::warning(this, tr("Cannot use this plate solution"), QString::fromStdString(e.what()));
        return false;
    }
}

// ------------------------------------------------------------------ export
void MainWindow::exportCandidates() {
    if (!m_result) return;
    const QString path = QFileDialog::getSaveFileName(this, tr("Export candidates"), "candidates.csv", tr("CSV (*.csv)"));
    if (path.isEmpty()) return;
    try {
        write_candidates_csv(path.toStdString(), *m_result);
        m_status->setText(tr("Wrote %1").arg(path));
    } catch (const std::exception& e) {
        QMessageBox::warning(this, tr("Export failed"), QString::fromStdString(e.what()));
    }
}

void MainWindow::exportTracklets() {
    if (!m_result || !m_prep) return;
    const QString path = QFileDialog::getSaveFileName(this, tr("Export tracklets"), "tracklets.csv", tr("CSV (*.csv)"));
    if (path.isEmpty()) return;
    try {
        write_tracklets_csv(path.toStdString(), *m_prep, *m_result, true);
        m_status->setText(tr("Wrote %1").arg(path));
    } catch (const std::exception& e) {
        QMessageBox::warning(this, tr("Export failed"), QString::fromStdString(e.what()));
    }
}

void MainWindow::about() {
    QMessageBox::about(this, tr("About FaintTrail"),
                       tr("<b>FaintTrail %1</b><br>Synthetic tracking engine for faint asteroids and comets.<br><br>"
                          "Shift-and-stack search over a grid of linear motions with a threshold calibrated on "
                          "time-shuffled copies of the data.<br><br>Free software under the GNU General Public License, version 3 or (at your option) any later version. There is NO WARRANTY.")
                           .arg(version()));
}
