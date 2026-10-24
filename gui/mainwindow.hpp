// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#pragma once

#include <QFutureWatcher>
#include <QMainWindow>
#include <QMutex>
#include <atomic>
#include <memory>
#include <vector>

#include "fainttrail/fainttrail.hpp"

class ImageView;
class QAction;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QProgressBar;
class QPushButton;
class QSlider;
class QSpinBox;
class QTableWidget;
class QTimer;

struct SearchOutcome {
    std::shared_ptr<fainttrail::Prepared> prep;
    std::shared_ptr<fainttrail::SearchResult> result;
    QString error;
    bool cancelled = false;
};

struct SimOutcome {
    std::shared_ptr<fainttrail::SimResult> sim;
    fainttrail::SimParams params;
    QString error;
};

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);

    // Load frames from FITS files (sorted by time).
    bool loadFiles(const QStringList& files);
    // Generate a synthetic sequence (asynchronously); `then_search` runs a search afterwards.
    void simulate(const fainttrail::SimParams& params, bool then_search);
    void runSearch();
    // For screenshots: called when an asynchronous search has finished.
    void setOnSearchFinished(std::function<void()> fn) { m_onSearchFinished = std::move(fn); }
    void selectCandidate(int row);
    bool importPlateSolutionFrom(const QString& path);

private slots:
    void openFiles();
    void openFolder();
    void simulateDialog();
    void cancelSearch();
    void exportCandidates();
    void importPlateSolution();
    void exportTracklets();
    void searchFinished();
    void simulationFinished();
    void pollProgress();
    void layerChanged();
    void frameChanged(int index);
    void toggleBlink();
    void blinkStep();
    void tableSelectionChanged();
    void about();

private:
    void buildUi();
    void setFrames(std::vector<fainttrail::Frame> frames, const QString& origin);
    void clearResults();
    void refreshImage();
    void refreshMarkers();
    void refreshTable();
    void refreshDetails();
    void updateRateLabels();
    void setBusy(bool busy, const QString& what);
    fainttrail::SearchParams currentParams() const;
    double pixelScale() const;
    int frameCount() const { return static_cast<int>(m_frames.size()); }

    // data
    std::vector<fainttrail::Frame> m_frames;
    std::vector<fainttrail::SimTruth> m_truth;     // known movers when simulated
    std::shared_ptr<fainttrail::Prepared> m_prep;
    std::shared_ptr<fainttrail::SearchResult> m_result;
    int m_selected = -1;

    // async work
    QFutureWatcher<SearchOutcome> m_searchWatcher;
    QFutureWatcher<SimOutcome> m_simWatcher;
    bool m_searchAfterSim = false;
    std::atomic<bool> m_cancel{false};
    QMutex m_progressMutex;
    double m_progress = 0.0;
    QString m_stage;
    QTimer* m_progressTimer = nullptr;
    std::function<void()> m_onSearchFinished;

    // widgets
    ImageView* m_view = nullptr;
    QComboBox* m_layer = nullptr;
    QSlider* m_frameSlider = nullptr;
    QLabel* m_frameLabel = nullptr;
    QPushButton* m_blinkButton = nullptr;
    QTimer* m_blinkTimer = nullptr;
    QTableWidget* m_table = nullptr;
    QLabel* m_cut[3] = {nullptr, nullptr, nullptr};
    QLabel* m_details = nullptr;
    QLabel* m_summary = nullptr;
    QLabel* m_status = nullptr;
    QLabel* m_hover = nullptr;
    QProgressBar* m_progressBar = nullptr;

    QDoubleSpinBox* m_vmin = nullptr;
    QDoubleSpinBox* m_vmax = nullptr;
    QLabel* m_rateLabel = nullptr;
    QDoubleSpinBox* m_fwhm = nullptr;
    QDoubleSpinBox* m_falseRate = nullptr;
    QSpinBox* m_nullRuns = nullptr;
    QDoubleSpinBox* m_threshold = nullptr;
    QSpinBox* m_threads = nullptr;
    QComboBox* m_backend = nullptr;
    QPushButton* m_searchButton = nullptr;
    QPushButton* m_cancelButton = nullptr;

    QAction* m_actExportCand = nullptr;
    QAction* m_actImportWcs = nullptr;
    QAction* m_actExportTrack = nullptr;
};
