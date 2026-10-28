// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// fainttrail-gui: minimal Qt 6 front end for the fainttrail library.
//
//   fainttrail-gui [FITS files or a folder]
//   fainttrail-gui --demo                 simulate a sequence and search it
//   fainttrail-gui --demo --screenshot F  ... then save a window screenshot and quit
#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QStyleFactory>
#include <QTimer>

#include "mainwindow.hpp"

static int run(QApplication& app);

int main(int argc, char** argv) {
    fainttrail::tune_allocator_for_images();
    QApplication app(argc, argv);
    // a windowed program has no console: show fatal errors instead of vanishing
    try {
        return run(app);
    } catch (const std::exception& e) {
        QMessageBox::critical(nullptr, "fainttrail-gui", QString("Fatal error: %1").arg(QString::fromLocal8Bit(e.what())));
    } catch (...) {
        QMessageBox::critical(nullptr, "fainttrail-gui", "Fatal error (unknown exception)");
    }
    return 1;
}

static int run(QApplication& app) {
    QApplication::setApplicationName("fainttrail-gui");
    QApplication::setApplicationVersion(fainttrail::version());
    QApplication::setStyle(QStyleFactory::create("Fusion"));

    QCommandLineParser cli;
    cli.setApplicationDescription("Synthetic tracking for faint asteroids and comets");
    cli.addHelpOption();
    cli.addVersionOption();
    QCommandLineOption demo("demo", "Simulate a sequence and search it.");
    QCommandLineOption shot("screenshot", "Save a screenshot after the search and quit.", "file");
    QCommandLineOption seed("seed", "Simulation seed for --demo.", "n", "1");
    QCommandLineOption wcsFrom("wcs-from", "Borrow the plate solution of a solved image.", "file");
    cli.addOption(wcsFrom);
    cli.addOption(demo);
    cli.addOption(shot);
    cli.addOption(seed);
    cli.addPositionalArgument("frames", "Registered FITS frames or a folder.", "[frames...]");
    cli.process(app);

    MainWindow w;
    w.show();

    QStringList files;
    for (const QString& a : cli.positionalArguments()) {
        QFileInfo fi(a);
        if (fi.isDir()) {
            for (const auto& e : QDir(a).entryInfoList({"*.fits", "*.fit", "*.fts"}, QDir::Files, QDir::Name))
                files << e.absoluteFilePath();
        } else {
            files << a;
        }
    }

    if (cli.isSet(shot)) {
        const QString path = cli.value(shot);
        w.setOnSearchFinished([&w, path]() {
            QTimer::singleShot(400, &w, [&w, path]() {
                w.grab().save(path);
                QApplication::quit();
            });
        });
    }

    if (!files.isEmpty()) {
        bool ok = w.loadFiles(files);
        if (ok && cli.isSet(wcsFrom)) ok = w.importPlateSolutionFrom(cli.value(wcsFrom));
        if (ok && cli.isSet(shot)) w.runSearch();
    } else if (cli.isSet(demo)) {
        fainttrail::SimParams p;
        p.width = 512;
        p.height = 384;
        p.seed = cli.value(seed).toULongLong();
        w.simulate(p, true);
    }
    return app.exec();
}
