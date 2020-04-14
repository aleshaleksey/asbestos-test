use crate::spectrum::Spectrum;
use crate::structs::XYDatum;

use std;
// use std::result::Result as R;
use std::io::BufRead;
// use std::ops::Try;

const CSV: &str = ".csv";

pub fn read_csv_to_spectrum(name: &str) -> std::io::Result<Spectrum> {
    // Check if the file is pretending to be a csv.
    if !name.to_lowercase().ends_with(CSV) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            "File does not end in csv. Will not read data.",
        ));
    }

    // read file to buffer, or at least plug it into one.
    let f = std::fs::File::open(name)?;
    let f = std::io::BufReader::new(f);

    let mut output_spectrum = Spectrum::new();

    // Read lines of CSV into the reader and make a spectrum.
    for line in f.lines().map(|l| l.expect("Expected string, got error.")) {
        // Turn the text line into a Vec<f64>
        let maybe_data: Vec<f64> = line
            .split(',')
            .map(|text| {
                text.parse::<f64>()
                    .expect("Could not parse csv input to numerical value")
            })
            .collect::<Vec<_>>();

        let l = maybe_data.len();
        if (l != 0) && (l != 2) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "File contains wrong number of columns (must be 2).",
            ));
        }

        if l == 0 {
            continue;
        }

        let d = XYDatum::new(maybe_data[0], maybe_data[1]);
        output_spectrum.add_point(d);
    }

    Ok(output_spectrum)
}
