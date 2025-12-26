// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef IPTSD_CONTACTS_DETECTION_DETECTOR_HPP
#define IPTSD_CONTACTS_DETECTION_DETECTOR_HPP

#include "../contact.hpp"
#include "algorithms/cluster.hpp"
#include "algorithms/convolution.hpp"
#include "algorithms/gaussian.hpp"
#include "algorithms/kernels.hpp"
#include "algorithms/maximas.hpp"
#include "algorithms/neutral.hpp"
#include "algorithms/overlaps.hpp"
#include "algorithms/watershed.hpp"
#include "config.hpp"

#include <common/casts.hpp>
#include <common/types.hpp>

#include <gsl/gsl>
#include <spdlog/spdlog.h>

#include <cmath>
#include <type_traits>
#include <vector>

namespace iptsd::contacts::detection {

template <class T, class TFit = T>
class Detector {
public:
	static_assert(std::is_floating_point_v<T>);
	static_assert(std::is_floating_point_v<TFit>);

private:
	Config<T> m_config;

	// The diagonal of the heatmap.
	T m_input_diagonal = casts::to<T>(0);

	// The heatmap with the neutral value subtracted.
	Image<T> m_img_neutral {};

	// Heatmap used for maxima + clustering + watershed (unblurred neutral-subtracted heatmap).
	Image<T> m_img_cluster {};

	// The blurred heatmap (used only for gaussian fitting).
	Image<T> m_img_blurred {};

	// The kernel that is used for blurring.
	Matrix3<T> m_kernel_blur = kernels::gaussian<T, 3, 3>(gsl::narrow_cast<T>(0.75));

	// The list of local maximas.
	std::vector<Point> m_maximas {};

	// The list of spanned clusters.
	std::vector<Box> m_clusters {};

	// Temporary storage for cluster merging.
	std::vector<Box> m_clusters_temp {};

	// Input parameters for gaussian fitting.
	std::vector<gaussian::Parameters<TFit>> m_fitting_params {};

	// Temporary storage for gaussian fitting.
	Image<TFit> m_fitting_temp {};

	// How many frames are left before the neutral value has to be recalculated.
	usize m_counter = 0;

	// The cached neutral value of the heatmap.
	T m_neutral = casts::to<T>(0);

public:
	Detector(Config<T> config) : m_config {std::move(config)} {};

	/*!
	 * Search for contacts in a capacitive heatmap.
	 *
	 * This function uses a threshold based recursive descent to build
	 * a list of connected clusters based on the (pre-processed) heatmap,
	 * and then uses gaussian fitting to fit an ellipse onto these clusters.
	 *
	 * @param[in] heatmap The heatmap to process.
	 * @param[out] contacts The list of detected contacts.
	 */
	template <int Rows, int Cols>
	void detect(const ImageBase<T, Rows, Cols> &heatmap, std::vector<Contact<T>> &contacts)
	{
		const Vector2<Eigen::Index> one = Vector2<Eigen::Index>::Ones();

		const Eigen::Index cols = heatmap.cols();
		const Eigen::Index rows = heatmap.rows();

		const Vector2<Eigen::Index> dimensions {cols - 1, rows - 1};

		const Eigen::Index bcols = m_img_neutral.cols();
		const Eigen::Index brows = m_img_neutral.rows();

		// Resize the internal buffers if neccessary.
		if (brows != rows || bcols != cols) {
			m_img_neutral.conservativeResize(rows, cols);
			m_img_cluster.conservativeResize(rows, cols);
			m_img_blurred.conservativeResize(rows, cols);
			m_fitting_temp.conservativeResize(rows, cols);

			if (m_config.normalize)
				m_input_diagonal = std::hypot(dimensions.x(), dimensions.y());
		}

		contacts.clear();
		m_clusters.clear();
		m_fitting_params.clear();
		m_maximas.clear();

		// Recalculate the neutral value if neccessary
		if (m_counter == 0) {
			m_neutral = neutral::calculate(heatmap,
			                               m_config.neutral_value_algorithm,
			                               m_config.neutral_value_offset);
		}

		// Update counter
		m_counter = (m_counter + 1) % m_config.neutral_value_backoff;

		// Subtract the neutral value from the whole heatmap
		m_img_neutral = (heatmap - m_neutral).max(casts::to<T>(0));

		// Unblurred clustering image
		m_img_cluster = m_img_neutral;

		// Blur only for gaussian fitting
		convolution::run(m_img_neutral, m_kernel_blur, m_img_blurred);

		const T athresh = m_config.activation_threshold;
		const T dthresh = m_config.deactivation_threshold;

		// Search for local maximas on the same image we use for clustering/splitting
		maximas::find(m_img_cluster, athresh, m_maximas);

		// Build initial clusters from the same image (avoid blur bridges)
		for (const Point &point : m_maximas) {
			Box cluster = cluster::span(m_img_cluster, point, athresh, dthresh);

			if (cluster.isEmpty())
				continue;

			// Extend the sides of the cluster by one pixel
			cluster.min() = (cluster.min() - one).cwiseMax(0);
			cluster.max() = (cluster.max() + one).cwiseMin(dimensions);

			// min() and max() are inclusive so we need to add one
			const Vector2<Eigen::Index> size = cluster.sizes() + one;

			// For gaussian fitting, the clusters should have at least 3x3 pixels
			if (size.x() < 3 || size.y() < 3)
				continue;

			m_clusters.push_back(std::move(cluster));
		}

		// Merge overlapping clusters
		overlaps::merge(m_clusters, m_clusters_temp, 5);
		
		std::vector<Box> final_clusters;

		// Watershed split
		if (m_config.detection_algorithm == Algorithm::WATERSHED) {
			spdlog::info("Using watershed algorithm. {} keep ratio, {} palm diagonal fraction.", m_config.watershed_keep_ratio, m_config.watershed_palm_diag_frac);
			final_clusters = {};
			final_clusters.reserve(m_clusters.size());

			// Only split if the blob has a small number of strong peaks (two-finger case)
			const usize max_ws_seeds = 3;

			// Palm gate: if bbox diagonal is large, do NOT split (preserve palm heuristics)
			// Tune 0.20..0.35. Higher => fewer splits on big blobs.
			const T palm_diag_frac = casts::to<T>(m_config.watershed_palm_diag_frac);
			const T heatmap_diag = std::hypot(casts::to<T>(dimensions.x()), casts::to<T>(dimensions.y()));
			const T palm_diag_px = heatmap_diag * palm_diag_frac;

			auto cluster_diag_px = [&](const Box &b) -> T {
				const auto sz = b.sizes() + one; // inclusive bbox size
				return std::hypot(casts::to<T>(sz.x()), casts::to<T>(sz.y()));
			};

			for (const Box &cluster : m_clusters) {
				// Preserve large blobs as-is (palm-ish)
				if (cluster_diag_px(cluster) >= palm_diag_px) {
					final_clusters.push_back(cluster);
					continue;
				}

				// Collect maxima inside merged cluster bbox
				std::vector<Point> seeds {};
				seeds.reserve(6);

				T peak_max = casts::to<T>(0);

				for (const Point &p : m_maximas) {
					if (!cluster.contains(p))
						continue;

					const T v = m_img_cluster(p.y(), p.x());
					if (v <= athresh)
						continue;

					peak_max = std::max(peak_max, v);
					seeds.push_back(p);
				}

				if (seeds.size() <= 1) {
					final_clusters.push_back(cluster);
					continue;
				}

				// Keep only strong seeds relative to the strongest peak
				if (peak_max > casts::to<T>(0)) {
					std::vector<Point> strong {};
					strong.reserve(seeds.size());

					for (const Point &p : seeds) {
						const T v = m_img_cluster(p.y(), p.x());
						if (v >= peak_max * m_config.watershed_keep_ratio)
							strong.push_back(p);
					}

					seeds.swap(strong);
				}

				// Too many strong peaks => likely palm/smear, preserve blob
				if (seeds.size() <= 1 || seeds.size() > max_ws_seeds) {
					final_clusters.push_back(cluster);
					continue;
				}

				auto splits = watershed::split(m_img_cluster, cluster, seeds, dthresh, 3);

				if (splits.size() < 2) {
					final_clusters.push_back(cluster);
					continue;
				}

				for (auto &b : splits)
					final_clusters.push_back(std::move(b));
			}
		} else {
			spdlog::info("Skipping watershed algorithm.");
			final_clusters = m_clusters;
		}

		// Prepare clusters for gaussian fitting.
		for (const Box &cluster : final_clusters) {
			const Vector2<TFit> mean = cluster.cast<TFit>().center();
			const Matrix2<TFit> prec = Matrix2<TFit>::Identity();

			const Vector2<Eigen::Index> size = cluster.sizes() + one;

			gaussian::Parameters<TFit> params {
				true,
				1,
				mean,
				prec,
				cluster,
				Image<TFit> {size.y(), size.x()},
			};

			m_fitting_params.push_back(std::move(params));
		}

		// Run gaussian fitting
		gaussian::fit(m_fitting_params, m_img_blurred, m_fitting_temp, 6);

		// Create a contact from every gaussian fitting parameter
		for (const auto &p : m_fitting_params) {
			if (!p.valid)
				continue;

			const Matrix2<TFit> cov = p.prec.inverse();

			Eigen::SelfAdjointEigenSolver<Matrix2<TFit>> solver {};
			solver.computeDirect(cov);

			Vector2<TFit> mean = p.mean;

			Vector2<TFit> size {};
			size.x() = std::sqrt(solver.eigenvalues().x()) * casts::to<TFit>(2);
			size.y() = std::sqrt(solver.eigenvalues().y()) * casts::to<TFit>(2);

			TFit orientation = std::atan2(solver.eigenvectors()(1, 1), solver.eigenvectors()(0, 1));
			
			// Normalize dimensions.
			if (m_config.normalize) {
				mean = (mean.array() / dimensions.cast<TFit>().array()).matrix();
				size = (size.array() / m_input_diagonal).matrix();
				orientation /= gsl::narrow_cast<TFit>(M_PI);
			}

			contacts.push_back(Contact<T> {mean.template cast<T>(),
			                               size.template cast<T>(),
			                               gsl::narrow_cast<T>(orientation),
			                               m_config.normalize});
		}
	}
};

} // namespace iptsd::contacts::detection

#endif // IPTSD_CONTACTS_DETECTION_DETECTOR_HPP
