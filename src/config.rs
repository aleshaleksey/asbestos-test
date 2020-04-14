use crate::writer::*;
use crate::{ACTINOLITE, AMOSITE, ANTHOPHYLLITE, CHRYSOTILE, CROCIDOLITE, TREMOLITE};
use serde;
use std;
use std::io::Read;
use toml;

#[derive(serde::Deserialize, serde::Serialize, Debug, Clone)]
pub struct AsbestosConfig {
    pub clean_source_dir: std::path::PathBuf,
    pub sample_source_dir: std::path::PathBuf,
    pub output_dir: std::path::PathBuf,
    pub amosite_file: Option<String>,
    pub crocidolite_file: Option<String>,
    pub chrysotile_file: Option<String>,
    pub actinolite_file: Option<String>,
    pub anthophyllite_file: Option<String>,
    pub tremolite_file: Option<String>,
    // The minimum width of the crocidolite-amosite peak.
    pub ca_min_width: f64,
    // The noise level of the machine in Abs.
    pub noise_level: f64,
}

impl AsbestosConfig {
    /// Create an asbestos config object from a file.
    pub fn from_file(path: &str) -> std::io::Result<AsbestosConfig> {
        let mut file = std::fs::File::open(path)?;
        let mut config_string = String::with_capacity(1000);
        file.read_to_string(&mut config_string)?;
        match toml::from_str(&config_string) {
            Ok(c) => Ok(c),
            Err(_) => Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "Config file could not be interpreted. Please check it for errors.",
            )),
        }
    }

    /// A function to run a function for each reference spectrum supplied
    /// in the config file.
    pub fn do_for_each_asbestos(&self, do_for: &dyn Fn(&Self, &mut Vec<AsbTypes>, &str, &str)) {
        let mut analysis_no = 0;
        let mut results = Vec::with_capacity(1000);
        if let Some(ref am) = self.amosite_file {
            do_for(&self, &mut results, &am, AMOSITE);
            analysis_no += 1;
        }
        if let Some(ref cr) = self.crocidolite_file {
            do_for(&self, &mut results, &cr, CROCIDOLITE);
            analysis_no += 1;
        }
        if let Some(ref ch) = self.chrysotile_file {
            do_for(&self, &mut results, &ch, CHRYSOTILE);
            analysis_no += 1;
        }
        if let Some(ref ac) = self.actinolite_file {
            do_for(&self, &mut results, &ac, ACTINOLITE);
            analysis_no += 1;
        }
        if let Some(ref an) = self.anthophyllite_file {
            do_for(&self, &mut results, &an, ANTHOPHYLLITE);
            analysis_no += 1;
        }
        if let Some(ref tr) = self.tremolite_file {
            do_for(&self, &mut results, &tr, TREMOLITE);
            analysis_no += 1;
        }

        if analysis_no == 0 {
            panic!("No asbestos reference files provided. Analysis cannot be completed.");
        }
        results.sort_by(|r, re| {
            r.sample
                .partial_cmp(&re.sample)
                .expect("Could not order results when writing file.")
        });
        match crate::writer::write_analysis_to_csv(&results, &self) {
            Ok(_) => println!("Analysis complete."),
            Err(e) => panic!("Could not write analysis file: {:?}", e),
        };
    }
}

#[test]
// NB this is not a proper automated test, it's just to have a gander.
fn test_config() {
    let path = "C:\\Users\\Vladimir Zholobenko\\alek-code\\asbestos-test\\config.toml";
    let mut file = std::fs::File::open(path).expect("Couldn't open config.");
    let mut config_string = String::new();
    file.read_to_string(&mut config_string)
        .expect("Couldn't read to string");
    println!("Config string:\n{}", config_string);
    let asbestos_config: AsbestosConfig =
        toml::from_str(&config_string).expect("We couldn't toml it.");

    println!("Config toml:\n{:#?}", asbestos_config);
}

#[test]
fn test_config_backwards() {
    let ac = AsbestosConfig {
        clean_source_dir: std::path::PathBuf::from("E:\\"),
        sample_source_dir: std::path::PathBuf::from("E:\\"),
        output_dir: std::path::PathBuf::from("E:\\"),
        amosite_file: Some("amo".into()),
        crocidolite_file: Some("croc".into()),
        chrysotile_file: Some("chrys".into()),
        actinolite_file: Some("acti".into()),
        anthophyllite_file: Some("anth".into()),
        tremolite_file: Some("trem".into()),
    };
    println!("This is ac: {:#?}", ac);

    let ac_str = toml::to_string(&ac).expect("couldn't untoml.");
    println!("This is ac  as string: \n{}", ac_str);
}
