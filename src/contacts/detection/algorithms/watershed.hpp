// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef IPTSD_CONTACTS_DETECTION_ALGORITHMS_WATERSHED_HPP
#define IPTSD_CONTACTS_DETECTION_ALGORITHMS_WATERSHED_HPP

#include <common/types.hpp>

#include <queue>
#include <vector>

namespace iptsd::contacts::detection::watershed {

namespace impl {

template <class T>
struct Node {
	T v;
	int label;
	Eigen::Index x;
	Eigen::Index y;
};

// Max-heap by value (process highest values first)
template <class T>
struct NodeCmp {
	bool operator()(const Node<T> &a, const Node<T> &b) const noexcept
	{
		return a.v < b.v;
	}
};

inline bool in_box(const Box &b, Eigen::Index x, Eigen::Index y) noexcept
{
	return x >= b.min().x() && x <= b.max().x() && y >= b.min().y() && y <= b.max().y();
}

} // namespace impl

/*!
 * Split a cluster bbox into per-seed sub-clusters using marker-based watershed.
 *
 * - Operates only inside cluster bounds.
 * - Considers only pixels with value > deactivation_threshold.
 * - Seeds are maxima points (must lie inside cluster and above activation threshold in caller).
 *
 * Returns 0..N boxes, one per surviving seed region.
 */
template <class Derived>
std::vector<Box> split(const DenseBase<Derived> &heatmap,
                       const Box &cluster,
                       const std::vector<Point> &seeds,
                       const typename DenseBase<Derived>::Scalar deactivation_threshold,
                       const Eigen::Index min_side_px = 3)
{
	using T = typename DenseBase<Derived>::Scalar;

	std::vector<Box> out {};
	out.reserve(seeds.size());

	if (cluster.isEmpty() || seeds.empty())
		return out;

	const Eigen::Index w = cluster.sizes().x() + 1; // inclusive bbox
	const Eigen::Index h = cluster.sizes().y() + 1;

	// label image local to bbox
	Image<int> labels {h, w};
	labels.setConstant(-1);

	auto lx = [&](Eigen::Index x) { return x - cluster.min().x(); };
	auto ly = [&](Eigen::Index y) { return y - cluster.min().y(); };

	std::priority_queue<impl::Node<T>, std::vector<impl::Node<T>>, impl::NodeCmp<T>> pq {};

	// Initialize with seeds
	for (int i = 0; i < (int)seeds.size(); i++) {
		const Eigen::Index sx = seeds[(usize)i].x();
		const Eigen::Index sy = seeds[(usize)i].y();

		if (!impl::in_box(cluster, sx, sy))
			continue;

		const T sv = heatmap(sy, sx);
		if (sv <= deactivation_threshold)
			continue;

		const Eigen::Index ix = lx(sx);
		const Eigen::Index iy = ly(sy);

		if (labels(iy, ix) != -1)
			continue; // two seeds on same pixel: ignore the later one

		labels(iy, ix) = i;
		pq.push({sv, i, sx, sy});
	}

	if (pq.empty())
		return out;

	// Grow regions from highest to lowest values
	while (!pq.empty()) {
		const auto cur = pq.top();
		pq.pop();

		const Eigen::Index x = cur.x;
		const Eigen::Index y = cur.y;

		// 4-neighborhood
		const Eigen::Index nx[4] = {x + 1, x - 1, x + 0, x + 0};
		const Eigen::Index ny[4] = {y + 0, y + 0, y + 1, y - 1};

		for (int k = 0; k < 4; k++) {
			const Eigen::Index xx = nx[k];
			const Eigen::Index yy = ny[k];

			if (!impl::in_box(cluster, xx, yy))
				continue;

			const T v = heatmap(yy, xx);
			if (v <= deactivation_threshold)
				continue;

			const Eigen::Index ix = lx(xx);
			const Eigen::Index iy = ly(yy);

			int &dst = labels(iy, ix);
			if (dst == -1) {
				dst = cur.label;
				pq.push({v, cur.label, xx, yy});
			} else if (dst != cur.label) {
				// Meeting front: leave existing label as-is.
				// This implicitly creates a boundary at the first-touch frontier.
				// (If you want an explicit "watershed line", set dst = -2 here and
				// skip later bbox accumulation of -2.)
			}
		}
	}

	// Build bbox per label
	out.assign((usize)seeds.size(), Box {});
	for (auto &b : out)
		b.setEmpty();

	for (Eigen::Index yy = cluster.min().y(); yy <= cluster.max().y(); yy++) {
		for (Eigen::Index xx = cluster.min().x(); xx <= cluster.max().x(); xx++) {
			const Eigen::Index ix = lx(xx);
			const Eigen::Index iy = ly(yy);

			const int lab = labels(iy, ix);
			if (lab < 0)
				continue;

			Box &b = out[(usize)lab];
			if (b.isEmpty()) {
				b.min() = Point {xx, yy};
				b.max() = Point {xx, yy};
			} else {
				b.extend(Point {xx, yy});
			}
		}
	}

	// Filter empties / too-small regions
	std::vector<Box> filtered {};
	filtered.reserve(out.size());
	for (auto &b : out) {
		if (b.isEmpty())
			continue;

		const auto sz = b.sizes() + Vector2<Eigen::Index>::Ones(); // inclusive
		if (sz.x() < min_side_px || sz.y() < min_side_px)
			continue;

		filtered.push_back(std::move(b));
	}

	return filtered;
}

} // namespace iptsd::contacts::detection::watershed

#endif // IPTSD_CONTACTS_DETECTION_ALGORITHMS_WATERSHED_HPP
