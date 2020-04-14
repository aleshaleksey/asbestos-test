#![feature(try_trait)]
#![allow(dead_code)]

extern crate serde;
extern crate toml;

mod baseline;
mod config;
mod reader;
mod spectrum;
mod structs;
mod writer;

const MIN_W: f64 = 4.0;
const MAX_W: f64 = 30.0;
const MIN_H: f64 = 0.0;
// Prominence of a peak as compared to the MAX-MIN of the spectrum.
const MIN_P: f64 = 0.01;
const DEFAULT_DIR: &str = "asbestos_config.toml";
// The tolerance for x-unit error when finding peaks.
const X_TOLERANCE: f64 = 3.0;

const AMOSITE: &str = "amosite";
const CROCIDOLITE: &str = "crocidolite";
const CHRYSOTILE: &str = "chrysotile";
const ACTINOLITE: &str = "actinolite";
const ANTHOPHYLLITE: &str = "anthophyllite";
const TREMOLITE: &str = "tremolite";

fn do_for_one(
    asb_conf: &config::AsbestosConfig,
    results: &mut Vec<writer::AsbTypes>,
    ref_name: &str,
    asb_type: &str,
) {
    let ref_path = asb_conf.clean_source_dir.clone().join(ref_name);
    let ref_path = ref_path.to_str().expect("Reference path cannot be read.");

    let spectral_noise = asb_conf.noise_level;
    let ca_hhw = asb_conf.ca_min_width;

    // read initial clean spectrum
    let mut clean_spectrum = match reader::read_csv_to_spectrum(ref_path) {
        Ok(s) => s,
        Err(e) => panic!("Error reading data to spectrum: {:?}", e),
    };

    // Process template spectrum.
    clean_spectrum.sort_by_x();
    clean_spectrum.subtract_minimum();
    clean_spectrum.flatten_baseline();
    clean_spectrum.find_peaks(MIN_W, MAX_W, MIN_H, MIN_P);
    clean_spectrum.normalise_peak_prominence();

    let mask = clean_spectrum.create_masked_from(&clean_spectrum);
    let mut minusmask = match clean_spectrum.subtract_spectrum(&mask) {
        Some(m) => m,
        None => panic!("spectral sub failed. No overlap or something?"),
    };
    minusmask.find_peaks(MIN_W, MAX_W, MIN_H, MIN_P);
    minusmask.normalise_peak_prominence();

    // Write out file.
    let mut out_path = asb_conf
        .output_dir
        .join(&asb_type)
        .to_str()
        .expect("Output path could not be converted.")
        .to_owned();
    out_path.push_str(asb_type);

    match writer::write_spectrum_to_csv(&out_path, "-mask-sub.csv", &minusmask) {
        Ok(_) => {
            println!("Output successfully written.");
        }
        Err(e) => {
            println!("Error writing file out: {:?}.", e);
        }
    }

    if !asb_conf.sample_source_dir.is_dir() {
        panic!(
            "The sample path ({:?}) is not a valid directory.",
            asb_conf.sample_source_dir
        );
    } else {
        let csv_l = std::ffi::OsString::from("csv");
        let csv_b = std::ffi::OsString::from("CSV");
        for entry in asb_conf
            .sample_source_dir
            .read_dir()
            .expect("Could not read dir.")
        {
            if let Ok(entry) = entry {
                let p = entry.path();
                let ext = p.extension();
                if entry.path().is_file() && (ext == Some(&csv_l) || ext == Some(&csv_b)) {
                    println!("We will process {:?}", entry.path());
                    let mut path = p.to_str().expect("path could not be converted.").to_owned();
                    path = path.replace("\\\\", "\\");
                    let mut spectrum = match reader::read_csv_to_spectrum(&path) {
                        Ok(s) => s,
                        Err(e) => panic!("Error reading data to spectrum: {:?}", e),
                    };

                    // Process dirty spectrum.
                    spectrum.sort_by_x();
                    //spectrum.subtract_minimum();
                    //spectrum.flatten_baseline();

                    //Create mask of dirty spectrum using clean:
                    let mask = spectrum.create_masked_from(&minusmask);

                    let mut subtracted = match spectrum.subtract_spectrum(&mask) {
                        Some(s) => s,
                        None => panic!("spectral sub failed. No overlap or something?"),
                    };

                    subtracted.find_peaks(MIN_W, MAX_W, MIN_H, MIN_P);
                    subtracted.normalise_peak_prominence();

                    // Write out file.
                    let name = p
                        .file_name()
                        .expect("The file should have a name")
                        .to_str()
                        .expect("path could not be converted")
                        .to_owned();
                    let mut out_path = asb_conf
                        .output_dir
                        .join(&name)
                        .to_str()
                        .expect("Output path could not be converted.")
                        .to_owned();
                    out_path.push_str(asb_type);

                    let (contains, extra) = match asb_type {
                        AMOSITE => minusmask.compare_peaks_ac_special(
                            &subtracted,
                            2,
                            3.0,
                            spectral_noise,
                            ca_hhw,
                        ),
                        ACTINOLITE => (
                            minusmask.compare_peaks_leniant(&subtracted, 2, 4.5, spectral_noise),
                            false,
                        ),
                        CHRYSOTILE => (
                            minusmask.compare_peaks_leniant(&subtracted, 1, 4.0, spectral_noise),
                            false,
                        ),
                        CROCIDOLITE => minusmask.compare_peaks_ac_special(
                            &subtracted,
                            2,
                            3.0,
                            spectral_noise,
                            ca_hhw,
                        ),
                        _ => (
                            minusmask.compare_peaks_leniant(&subtracted, 2, 3.0, spectral_noise),
                            false,
                        ),
                    };

                    // Prepare a result.
                    let mut res = if results.iter().any(|r| r.sample == name) {
                        results.iter_mut().find(|r| r.sample == name).unwrap()
                    } else {
                        results.push(writer::AsbTypes::with_name(&name));
                        results.last_mut().unwrap()
                    };

                    match asb_type {
                        AMOSITE => {
                            res.amosite = contains;
                            if contains {
                                res.crocidolite = extra;
                            }
                        }
                        ACTINOLITE => {
                            res.actinolite = contains;
                        }
                        CHRYSOTILE => {
                            res.chrysotile = contains;
                        }
                        CROCIDOLITE => {
                            res.crocidolite = contains;
                            if contains {
                                res.amosite = extra;
                            }
                        }
                        TREMOLITE => {
                            res.tremolite = contains;
                        }
                        ANTHOPHYLLITE => {
                            res.anthophyllite = contains;
                        }
                        _ => panic!("Unknown asbestos type!. Impossible!"),
                    };

                    // slightly better test output.
                    println!("Sample {}. Contains {}: {}.", name, asb_type, contains);

                    match writer::write_spectrum_to_csv(&out_path, "-mask.csv", &mask) {
                        Ok(_) => {
                            println!("Mask successfully written.");
                        }
                        Err(e) => {
                            println!("Error writing file out: {:?}.", e);
                        }
                    }
                    match writer::write_spectrum_to_csv(&out_path, "-mask-sub.csv", &subtracted) {
                        Ok(_) => {
                            println!("Unmasked peaks successfully written.");
                        }
                        Err(e) => {
                            println!("Error writing file out: {:?}.", e);
                        }
                    }
                }
            }
        }
    }
}

fn main() {
    println!("Hello, asbestos!");
    let args = std::env::args().collect::<Vec<_>>();
    if args.len() == 1 {
        let config = match config::AsbestosConfig::from_file(DEFAULT_DIR) {
            Ok(c) => c,
            Err(e) => panic!(
                "No arguments supplied and \"asbestos_config.toml\" could not be read: {:?}",
                e,
            ),
        };
        //let function = |x| {println!("This is x={}",x);};
        config.do_for_each_asbestos(&do_for_one);
    } else if args.len() == 2 {
        let file_name = &args.last().expect("There should be at least one argument.");

        // read initial spectrum.
        let mut spectrum = match reader::read_csv_to_spectrum(file_name) {
            Ok(s) => s,
            Err(e) => panic!("Error reading data to spectrum: {:?}", e),
        };

        // Process spectrum.
        spectrum.sort_by_x();
        let s2 = spectrum.clone();
        spectrum.subtract_minimum();
        spectrum.flatten_baseline();
        spectrum.find_peaks(MIN_W, MAX_W, MIN_H, MIN_P);
        spectrum.normalise_peak_prominence();

        let spectrum_sub = s2.subtract_spectrum(&spectrum).unwrap();

        //Open the UNKNOWN spectrum

        // View data onscreen.
        println!(
            "The spectral data, point by point, look like this:\n {:#?}",
            spectrum
        );

        // write processed data as a csv.
        match writer::write_spectrum_to_csv(file_name, "-processed.csv", &spectrum) {
            Ok(_) => {
                println!("Output successfully written.");
            }
            Err(e) => {
                println!("Error writing file out: {:?}.", e);
            }
        }

        // write linear sub as csv.
        match writer::write_spectrum_to_csv(file_name, "-linear sub.csv", &spectrum_sub) {
            Ok(_) => {
                println!("Output successfully written.");
            }
            Err(e) => {
                println!("Error writing file out: {:?}.", e);
            }
        }

        let mask = spectrum.create_masked_from(&spectrum);

        // write masked as csv.
        match writer::write_spectrum_to_csv(file_name, "-mask_self.csv", &mask) {
            Ok(_) => {
                println!("Output successfully written.");
            }
            Err(e) => {
                println!("Error writing file out: {:?}.", e);
            }
        }

        let minusmask = match spectrum.subtract_spectrum(&mask) {
            Some(m) => m,
            None => panic!("spectral sub failed. No overlap or something?"),
        };

        // write masked as csv.
        match writer::write_spectrum_to_csv(file_name, "-mask-sub.csv", &minusmask) {
            Ok(_) => {
                println!("Output successfully written.");
            }
            Err(e) => {
                println!("Error writing file out: {:?}.", e);
            }
        }
    } else if args.len() == 3 {
        let file_name = &args[1];
        let mut clean_spectrum = match reader::read_csv_to_spectrum(file_name) {
            Ok(s) => s,
            Err(e) => panic!("Error reading data to spectrum: {:?}", e),
        };

        // Process template spectrum.
        clean_spectrum.sort_by_x();
        clean_spectrum.subtract_minimum();
        clean_spectrum.flatten_baseline();
        clean_spectrum.find_peaks(MIN_W, MAX_W, MIN_H, MIN_P);
        clean_spectrum.normalise_peak_prominence();

        // read initial spectrum.
        let file_name2 = &args.last().expect("Should be there because we checked.");
        let mut spectrum = match reader::read_csv_to_spectrum(file_name2) {
            Ok(s) => s,
            Err(e) => panic!("Error reading data to spectrum: {:?}", e),
        };

        // Process dirty spectrum.
        spectrum.sort_by_x();

        //Create mask of dirty spectrum using clean:
        let mask = spectrum.create_masked_from(&clean_spectrum);

        let mut subtracted = match spectrum.subtract_spectrum(&mask) {
            Some(s) => s,
            None => panic!("spectral sub failed. No overlap or something?"),
        };

        spectrum.subtract_minimum();
        spectrum.flatten_baseline();
        subtracted.find_peaks(MIN_W, MAX_W, MIN_H, MIN_P);
        subtracted.normalise_peak_prominence();

        // Print on screen comparison.
        println!(
            "Peaks of {:?} spectrum: {:#?}",
            file_name, clean_spectrum.peaks
        );
        println!(
            "Peaks of {:?} spectrum: {:#?}",
            file_name2, subtracted.peaks
        );

        // Write out file.
        match writer::write_spectrum_to_csv(file_name2, "-mask.csv", &mask) {
            Ok(_) => {
                println!("Output successfully written.");
            }
            Err(e) => {
                println!("Error writing file out: {:?}.", e);
            }
        }
        match writer::write_spectrum_to_csv(file_name2, "-mask-sub.csv", &subtracted) {
            Ok(_) => {
                println!("Output successfully written.");
            }
            Err(e) => {
                println!("Error writing file out: {:?}.", e);
            }
        }
    } else {
        println!("Wrong number of arguments. Stop arguing.");
    }
}
