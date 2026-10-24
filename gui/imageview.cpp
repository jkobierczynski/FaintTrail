// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
#include "imageview.hpp"

#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>

QImage to_qimage(const fainttrail::Image& img, Stretch stretch, double lo, double hi) {
    if (img.empty()) return {};
    if (!(lo < hi)) {
        std::vector<float> s;
        const int step = std::max(1, static_cast<int>(std::sqrt(img.size() / 40000.0)));
        for (int y = 0; y < img.height; y += step)
            for (int x = 0; x < img.width; x += step) {
                const float v = img.at(x, y);
                if (std::isfinite(v)) s.push_back(v);
            }
        if (s.empty()) return {};
        std::sort(s.begin(), s.end());
        lo = s[static_cast<std::size_t>(0.25 * (s.size() - 1))];
        const double med = s[s.size() / 2];
        hi = s[static_cast<std::size_t>(0.998 * (s.size() - 1))];
        // put the sky a little above black
        lo = med - 2.0 * (med - lo);
        if (!(hi > lo)) hi = lo + 1.0;
    }
    QImage out(img.width, img.height, QImage::Format_Grayscale8);
    const double k = 12.0, norm = std::asinh(k);
    for (int y = 0; y < img.height; ++y) {
        uchar* row = out.scanLine(y);
        for (int x = 0; x < img.width; ++x) {
            const float v = img.at(x, y);
            double u = std::isfinite(v) ? (v - lo) / (hi - lo) : 0.0;
            u = std::clamp(u, 0.0, 1.0);
            if (stretch == Stretch::Asinh) u = std::asinh(k * u) / norm;
            row[x] = static_cast<uchar>(std::lround(255.0 * u));
        }
    }
    return out;
}

ImageView::ImageView(QWidget* parent) : QWidget(parent) {
    setMouseTracking(true);
    setMinimumSize(320, 240);
    setAutoFillBackground(false);
    setFocusPolicy(Qt::StrongFocus);
}

void ImageView::setImage(const fainttrail::Image& img, Stretch stretch, double lo, double hi) {
    const bool size_changed = img.width != m_data.width || img.height != m_data.height;
    m_data = img;
    m_qimage = to_qimage(img, stretch, lo, hi);
    if (size_changed || m_fitPending) fitToWindow();
    update();
}

void ImageView::setMarkers(std::vector<Marker> markers) {
    m_markers = std::move(markers);
    update();
}

void ImageView::fitToWindow() {
    if (m_data.empty() || width() <= 0) {
        m_fitPending = true;
        return;
    }
    m_fitPending = false;
    m_scale = std::min(static_cast<double>(width()) / m_data.width, static_cast<double>(height()) / m_data.height);
    m_offset = QPointF((width() - m_data.width * m_scale) / 2.0, (height() - m_data.height * m_scale) / 2.0);
    update();
}

void ImageView::centerOn(double x, double y) {
    if (m_data.empty()) return;
    m_scale = std::max(m_scale, 3.0);
    m_offset = QPointF(width() / 2.0 - (x + 0.5) * m_scale, height() / 2.0 - (y + 0.5) * m_scale);
    update();
}

QPointF ImageView::toImage(const QPointF& w) const { return (w - m_offset) / m_scale; }

QPointF ImageView::toWidget(double x, double y) const {
    // pixel centres are at integer coordinates; pixel (0,0) spans [0,1) in display space
    return m_offset + QPointF((x + 0.5) * m_scale, (y + 0.5) * m_scale);
}

void ImageView::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.fillRect(rect(), QColor(18, 20, 24));
    if (m_qimage.isNull()) {
        p.setPen(QColor(150, 155, 165));
        p.drawText(rect(), Qt::AlignCenter, tr("Open registered FITS frames or simulate a sequence"));
        return;
    }
    p.setRenderHint(QPainter::SmoothPixmapTransform, m_scale < 1.0);
    const QRectF target(m_offset.x(), m_offset.y(), m_data.width * m_scale, m_data.height * m_scale);
    p.drawImage(target, m_qimage);

    p.setRenderHint(QPainter::Antialiasing, true);
    for (const auto& m : m_markers) {
        const QPointF c = toWidget(m.x, m.y);
        QColor col = m.color;
        if (!m.track.empty() && m.selected) {
            QPainterPath path;
            path.moveTo(toWidget(m.track.front().x(), m.track.front().y()));
            for (const auto& t : m.track) path.lineTo(toWidget(t.x(), t.y()));
            QColor tc = col;
            tc.setAlpha(150);
            p.setPen(QPen(tc, 1.2, Qt::DashLine));
            p.drawPath(path);
        }
        const double r = std::max(5.0, m.radius * std::min(m_scale, 4.0));
        p.setPen(QPen(col, m.selected ? 2.4 : 1.4));
        p.setBrush(Qt::NoBrush);
        if (m.cross) {
            p.drawLine(c + QPointF(-r, 0), c + QPointF(-r / 3, 0));
            p.drawLine(c + QPointF(r / 3, 0), c + QPointF(r, 0));
            p.drawLine(c + QPointF(0, -r), c + QPointF(0, -r / 3));
            p.drawLine(c + QPointF(0, r / 3), c + QPointF(0, r));
        } else {
            p.drawEllipse(c, r, r);
        }
        if (!m.label.isEmpty()) {
            p.setPen(col);
            p.drawText(c + QPointF(r + 3, -r + 2), m.label);
        }
    }
}

void ImageView::wheelEvent(QWheelEvent* e) {
    if (m_data.empty()) return;
    const QPointF pos = e->position();
    const QPointF before = toImage(pos);
    const double factor = std::pow(1.0015, e->angleDelta().y());
    m_scale = std::clamp(m_scale * factor, 0.05, 40.0);
    m_offset = pos - before * m_scale;
    update();
}

void ImageView::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) {
        // click on a marker selects it
        for (int i = static_cast<int>(m_markers.size()) - 1; i >= 0; --i) {
            const auto& m = m_markers[static_cast<std::size_t>(i)];
            if (m.cross) continue;
            const QPointF c = toWidget(m.x, m.y);
            const double r = std::max(5.0, m.radius * std::min(m_scale, 4.0)) + 3.0;
            if (QLineF(c, e->position()).length() <= r) {
                emit markerClicked(i);
                return;
            }
        }
        m_dragging = true;
        m_dragStart = e->position();
        m_offsetStart = m_offset;
        setCursor(Qt::ClosedHandCursor);
    }
}

void ImageView::mouseMoveEvent(QMouseEvent* e) {
    if (m_dragging) {
        m_offset = m_offsetStart + (e->position() - m_dragStart);
        update();
    }
    if (!m_data.empty()) {
        const QPointF ip = toImage(e->position()) - QPointF(0.5, 0.5);
        const int x = static_cast<int>(std::lround(ip.x())), y = static_cast<int>(std::lround(ip.y()));
        const bool inside = m_data.inside(x, y);
        emit hovered(ip.x(), ip.y(), inside ? m_data.at(x, y) : 0.0, inside);
    }
}

void ImageView::mouseReleaseEvent(QMouseEvent*) {
    m_dragging = false;
    unsetCursor();
}

void ImageView::mouseDoubleClickEvent(QMouseEvent*) { fitToWindow(); }

void ImageView::resizeEvent(QResizeEvent*) {
    if (m_fitPending) fitToWindow();
}
