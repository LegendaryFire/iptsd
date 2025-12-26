// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef IPTSD_CONTACTS_DETECTION_CONFIG_HPP
#define IPTSD_CONTACTS_DETECTION_CONFIG_HPP

#include "algorithms/neutral.hpp"

#include <common/casts.hpp>
#include <common/types.hpp>

namespace iptsd::contacts::detection {
enum class Algorithm : u8 {
	LEGACY,
	WATERSHED,
};

template <class T>
struct Config {
public:
	static_assert(std::is_floating_point_v<T>);

public:
	/*
	 * Whether output dimensions should be normalized.
	 */
	bool normalize = false;

	/*
	 * How the neutral value of the heatmap is calculated.
	 */
	enum neutral::Algorithm neutral_value_algorithm = neutral::Algorithm::MODE;

	/*
	 * Configuration to enable or disable watershed.
	 */
	enum detection::Algorithm detection_algorithm = detection::Algorithm::LEGACY;

	/*
	 * An offset that is added to the calculated neutral value.
	 * If neutral_value_algorithm is set to CONSTANT, this defines the neutral value.
	 */
	T neutral_value_offset = casts::to<T>(0);

	/*
	 * How many frames to wait before recalculating the neutral value.
	 * A value of 1 means to recalculate the neutral value every frame.
	 */
	usize neutral_value_backoff = 1;

	/*
	 * If a pixel of the input data is larger than this value plus the neutral value
	 * it is marked as a contact and a recursive cluster search is started.
	 */
	T activation_threshold = casts::to<T>(24);

	/*
	 * If a pixel of the input data is below this value plus the neutral value,
	 * the recursive cluster search will stop once it reaches it.
	 */
	T deactivation_threshold = casts::to<T>(20);

	/*
	 * Configuration parameter to fine-tune watershed.
	 */
	T watershed_keep_ratio = casts::to<T>(0.70);

	/*
	 * Configuration parameter to fine-tune watershed.
	 */
	T watershed_palm_diag_frac = casts::to<T>(0.25);
};

} // namespace iptsd::contacts::detection

#endif // IPTSD_CONTACTS_DETECTION_CONFIG_HPP
