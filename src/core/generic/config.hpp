// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef IPTSD_CORE_GENERIC_CONFIG_HPP
#define IPTSD_CORE_GENERIC_CONFIG_HPP

#include "errors.hpp"

#include <common/error.hpp>
#include <common/types.hpp>
#include <contacts/config.hpp>
#include <ipts/parser.hpp>

#include <optional>
#include <string>

namespace iptsd::core {

class Config {
public:
	// [Config]
	bool invert_x = false;
	bool invert_y = false;

	f64 width = 0;
	f64 height = 0;

	// [Touchscreen]
	bool touchscreen_disable = false;
	bool touchscreen_disable_on_palm = false;
	bool touchscreen_disable_on_stylus = false;
	f64 touchscreen_overshoot = 0.5;

	// [Touchpad]
	bool touchpad_disable = false;
	bool touchpad_disable_on_palm = false;
	f64 touchpad_overshoot = 0.5;

	// [Contacts]
	std::string contacts_neutral = "mode";
	f64 contacts_neutral_value = 0;
	f64 contacts_activation_threshold = 40;
	f64 contacts_deactivation_threshold = 36;
	f64 contacts_size_thresh_min = 0.1;
	f64 contacts_size_thresh_max = 0.5;
	f64 contacts_position_thresh_min = 0.04;
	f64 contacts_position_thresh_max = 2;
	f64 contacts_orientation_thresh_min = 1;
	f64 contacts_orientation_thresh_max = 15;
	f64 contacts_size_min = 0.2;
	f64 contacts_size_max = 2;
	f64 contacts_aspect_min = 1;
	f64 contacts_aspect_max = 2.5;

	// [MultiTouch]
	f64 multitouch_neutral_value = 0;
	f64 multitouch_activation_threshold = 40;
	f64 multitouch_deactivation_threshold = 36;

	// [Stylus]
	bool stylus_disable = false;
	f64 stylus_tip_distance = 0;

	// [DFT]
	usize dft_position_min_amp = 50;
	usize dft_position_min_mag = 2000;
	f64 dft_position_exp = -0.7;
	usize dft_button_min_mag = 1000;
	usize dft_freq_min_mag = 10000;
	usize dft_tilt_min_mag = 10000;
	usize dft_mpp2_button_min_mag = 50000;
	usize dft_mpp2_contact_min_mag = 50000;
	f64 dft_tilt_distance = 0.6;
	bool dft_allow_split_events = false;

public:
	/*!
	 * Generates a configuration object for the contact detection library.
	 *
	 * @return A config object for contact detection using doubles.
	 */
	[[nodiscard]] contacts::Config<f64> contacts() const
	{
		contacts::Config<f64> config {};

		config.detection.normalize = true;
		config.detection.contacts_activation_threshold = this->contacts_activation_threshold / 255.0;
		config.detection.contacts_deactivation_threshold = this->contacts_deactivation_threshold / 255.0;
		config.detection.multitouch_activation_threshold = this->multitouch_activation_threshold / 255.0;
		config.detection.multitouch_deactivation_threshold = this->multitouch_deactivation_threshold / 255.0;

		using Algorithm = contacts::detection::neutral::Algorithm;

		if (this->contacts_neutral == "mode")
			config.detection.neutral_value_algorithm = Algorithm::MODE;
		else if (this->contacts_neutral == "average")
			config.detection.neutral_value_algorithm = Algorithm::AVERAGE;
		else if (this->contacts_neutral == "constant")
			config.detection.neutral_value_algorithm = Algorithm::CONSTANT;
		else
			throw common::Error<Error::InvalidNeutralValueAlgorithm> {};

		const f64 nval_contacts = this->contacts_neutral_value;
		const f64 nval_multitouch = this->multitouch_neutral_value;


		config.detection.contacts_neutral_offset = nval_contacts / 255.0;
		config.detection.multitouch_neutral_offset = nval_multitouch / 255.0;
		config.detection.neutral_value_backoff = 16; // TODO: config option

		const f64 diagonal = std::hypot(this->width, this->height);

		config.validation.track_validity = true;
		config.validation.size_limits = Vector2<f64> {
			this->contacts_size_min / diagonal,
			this->contacts_size_max / diagonal,
		};
		config.validation.aspect_limits = Vector2<f64> {
			this->contacts_aspect_min,
			this->contacts_aspect_max,
		};

		config.stability.size_threshold = Vector2<f64> {
			this->contacts_size_thresh_min / diagonal,
			this->contacts_size_thresh_max / diagonal,
		};
		config.stability.position_threshold = Vector2<f64> {
			this->contacts_position_thresh_min / diagonal,
			this->contacts_position_thresh_max / diagonal,
		};
		config.stability.orientation_threshold = Vector2<f64> {
			this->contacts_orientation_thresh_min / 180,
			this->contacts_orientation_thresh_max / 180,
		};

		return config;
	}
};

} // namespace iptsd::core

#endif // IPTSD_CORE_GENERIC_CONFIG_HPP
