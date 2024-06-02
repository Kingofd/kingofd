/*
* This is a simple trait generation algorithm that ensures that certain
* traits are only given to a unit when all required prerequisite traits
* are present in the unit or the other way around in the case of 
* excluding traits. This was written during the Open Source Lab sofware
* practical at TUM (Technical University of Munich). It can be found in 
* its entirety in the Battle for Wesnoth repository.
*/

void unit::generate_traits(bool must_have_only)
{
	// [...]
	std::vector<const config*> candidate_traits;
	std::vector<std::string> temp_require_traits;
	std::vector<std::string> temp_exclude_traits;

	// Now randomly fill out to the number of traits required or until
	// there aren't any more traits.
	int nb_traits = current_traits.size();
	int max_traits = u_type.num_traits();
	for(; nb_traits < max_traits; ++nb_traits)
	{
		current_traits = modifications_.child_range("trait");
		candidate_traits.clear();
		for(const config& t : u_type.possible_traits()) {
			// Skip the trait if the unit already has it.
			const std::string& tid = t["id"];
			bool already = false;
			for(const config& mod : current_traits) {
				if(mod["id"] == tid) {
					already = true;
					break;
				}
			}

			if(already) {
				continue;
			}
			// Skip trait if trait requirements are not met
			// or trait exclusions are present
			temp_require_traits = utils::split(t["require_traits"]);
			temp_exclude_traits = utils::split(t["exclude_traits"]);

			// See if the unit already has a trait that excludes the current one
			for(const config& mod : current_traits) {
				if (mod["exclude_traits"] != "") {
					for (const auto& c: utils::split(mod["exclude_traits"])) {
						temp_exclude_traits.push_back(c);
					}
				}
			}

			// First check for requirements
			bool trait_req_met = true;
			for(const std::string& s : temp_require_traits) {
				bool has_trait = false;
				for(const config& mod : current_traits) {
					if (mod["id"] == s)
						has_trait = true;
				}
				if(!has_trait) {
					trait_req_met = false;
					break;
				}
			}
			if(!trait_req_met) {
				continue;
			}

			// Now check for exclusionary traits
			bool trait_exc_met = true;

			for(const std::string& s : temp_exclude_traits) {
				bool has_exclusionary_trait = false;
				for(const config& mod : current_traits) {
					if (mod["id"] == s)
						has_exclusionary_trait = true;
				}
				if (tid == s) {
					has_exclusionary_trait = true;
				}
				if(has_exclusionary_trait) {
					trait_exc_met = false;
					break;
				}
			}
			if(!trait_exc_met) {
				continue;
			}

			const std::string& avl = t["availability"];
			// The trait is still available, mark it as a candidate for randomizing.
			// For leaders, only traits with availability "any" are considered.
			if(!must_have_only && (!can_recruit() || avl == "any")) {
				candidate_traits.push_back(&t);
			}
		}
		// No traits available anymore? Break
		if(candidate_traits.empty()) {
			break;
		}

		int num = randomness::generator->get_random_int(0,candidate_traits.size()-1);
		modifications_.add_child("trait", *candidate_traits[num]);
		candidate_traits.erase(candidate_traits.begin() + num);
	}
	// Once random traits are added, don't do it again.
	// Such as when restoring a saved character.
	random_traits_ = false;
}
