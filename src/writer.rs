use crate::config;
use crate::spectrum::Spectrum;
use crate::{ACTINOLITE, AMOSITE, ANTHOPHYLLITE, CHRYSOTILE, CROCIDOLITE, TREMOLITE};
use std;
use std::io::Write;

// A struct that holds the types of asbestos that were detected.
#[derive(Debug, Clone)]
pub struct AsbTypes {
    pub sample: String,
    pub amosite: bool,
    pub crocidolite: bool,
    pub chrysotile: bool,
    pub actinolite: bool,
    pub anthophyllite: bool,
    pub tremolite: bool,
}

impl AsbTypes {
    pub fn with_name(name: &str) -> AsbTypes {
        AsbTypes {
            sample: name.to_owned(),
            amosite: false,
            crocidolite: false,
            chrysotile: false,
            actinolite: false,
            anthophyllite: false,
            tremolite: false,
        }
    }
}

pub fn write_spectrum_to_csv(
    base_name: &str,
    suffix: &str,
    spectrum: &Spectrum,
) -> std::io::Result<()> {
    let mut name = base_name.to_owned().to_lowercase();
    name = name.replace(".csv", "");
    name.push_str(suffix);

    let file = std::fs::File::create(&name)?;
    let mut file = std::io::LineWriter::new(file);

    for datum in spectrum.data.iter() {
        let line = format!("{},{}\n", datum.x, datum.y);
        file.write_all(line.as_bytes())?;
    }
    file.flush()?;
    Ok(())
}

// A function to write the final analysis to a nice CSV file.
pub fn write_analysis_to_csv(
    results: &[AsbTypes],
    asb_c: &config::AsbestosConfig,
) -> std::io::Result<()> {
    let name = asb_c
        .output_dir
        .join("analysis.csv")
        .to_str()
        .expect("Cannot convert file name to real file name.")
        .to_string();
    let file = std::fs::File::create(&name)?;
    let mut file = std::io::LineWriter::new(file);

    let header = format!(
        "Sample,{},{},{},{},{},{}\n",
        ACTINOLITE, AMOSITE, ANTHOPHYLLITE, CHRYSOTILE, CROCIDOLITE, TREMOLITE,
    );
    file.write_all(header.as_bytes())?;

    for result in results.iter() {
        let act = if asb_c.actinolite_file.is_some() {
            result.actinolite.to_string()
        } else {
            "Not tested".to_string()
        };
        let amo = if asb_c.amosite_file.is_some() {
            result.amosite.to_string()
        } else {
            "Not tested".to_string()
        };
        let ant = if asb_c.anthophyllite_file.is_some() {
            result.anthophyllite.to_string()
        } else {
            "Not tested".to_string()
        };
        let chr = if asb_c.chrysotile_file.is_some() {
            result.chrysotile.to_string()
        } else {
            "Not tested".to_string()
        };
        let cro = if asb_c.crocidolite_file.is_some() {
            result.crocidolite.to_string()
        } else {
            "Not tested".to_string()
        };
        let tre = if asb_c.tremolite_file.is_some() {
            result.tremolite.to_string()
        } else {
            "Not tested".to_string()
        };

        let line = format!(
            "{},{},{},{},{},{},{}\n",
            result.sample, act, amo, ant, chr, cro, tre,
        );
        file.write_all(line.as_bytes())?;
    }
    file.flush()?;
    Ok(())
}
