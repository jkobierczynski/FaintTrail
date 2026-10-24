// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Zoomable, pannable image display with candidate overlays.
#pragma once

#include <QColor>
#include <QImage>
#include <QPointF>
#include <QString>
#include <QWidget>
#include <vector>

#include "fainttrail/types.hpp"

struct Marker {
    double x = 0, y = 0;        // 0-based image coordinates
    QColor color;
    double radius = 6;
    QString label;
    bool selected = false;
    bool cross = false;          // draw a cross instead of a circle
    std::vector<QPointF> track;  // optional predicted path
};

enum class Stretch { Asinh, Linear };

// Convert a float image into an 8-bit display image.
// lo/hi: data range; if lo >= hi a robust percentile range is used.
QImage to_qimage(const fainttrail::Image& img, Stretch stretch, double lo = 0, double hi = 0);

class ImageView : public QWidget {
    Q_OBJECT
public:
    explicit ImageView(QWidget* parent = nullptr);

    void setImage(const fainttrail::Image& img, Stretch stretch, double lo = 0, double hi = 0);
    void setMarkers(std::vector<Marker> markers);
    void centerOn(double x, double y);
    void fitToWindow();
    bool hasImage() const { return !m_qimage.isNull(); }

signals:
    void hovered(double x, double y, double value, bool inside);
    void markerClicked(int index);

protected:
    void paintEvent(QPaintEvent*) override;
    void wheelEvent(QWheelEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;
    void mouseDoubleClickEvent(QMouseEvent*) override;
    void resizeEvent(QResizeEvent*) override;

private:
    QPointF toImage(const QPointF& widget) const;
    QPointF toWidget(double x, double y) const;

    QImage m_qimage;
    fainttrail::Image m_data;
    std::vector<Marker> m_markers;
    double m_scale = 1.0;
    QPointF m_offset;            // widget position of image pixel (0,0) corner
    bool m_dragging = false;
    bool m_fitPending = true;
    QPointF m_dragStart, m_offsetStart;
};
