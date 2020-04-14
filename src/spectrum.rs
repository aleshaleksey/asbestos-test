use crate::structs::{MinMax, Peak, XYDatum};
use crate::X_TOLERANCE;

/// A structure that holds the data for a single spectrum.
/// NB, all x-axis data is stored for simplicity if segments of data need
/// to be inserted or removed.
#[derive(Debug, Clone)]
pub struct Spectrum {
    pub data: Vec<XYDatum>,
    pub peaks: Vec<Peak>,
}

impl Spectrum {
    /// Function to make a new instance of spectrum.
    pub fn new() -> Spectrum {
        Spectrum {
            data: Vec::with_capacity(100_000),
            peaks: Vec::new(),
        }
    }

    /// A function to add a data point into the spectrum.
    pub fn add_point(&mut self, d: XYDatum) {
        self.data.push(d)
    }

    /// A basic sort to organise points by x-value.
    pub fn sort_by_x(&mut self) {
        self.data.sort_by(|a, b| {
            a.x.partial_cmp(&b.x)
                .expect("Failed when sorting because x-data contained NaN.")
        });
    }

    /// Subtracts the minimum y height from the spectrum.
    pub fn subtract_minimum(&mut self) {
        let min = self.data.iter().fold(
            std::f64::MAX,
            |acc, datum| if datum.y < acc { datum.y } else { acc },
        );
        for d in self.data.iter_mut() {
            d.y -= min;
        }
    }

    /// A basic function to flatten a spectrum linearly vs edge points.
    pub fn flatten_baseline_inner(&mut self, range: &MinMax) {
        let l = self.data.len();
        let step_size = (range.max - range.min) / l as f64;

        for (i, dat) in self.data.iter_mut().enumerate() {
            dat.y -= step_size * i as f64;
        }
    }

    /// A function that flattens the baseline of a whole spectrum.
    pub fn flatten_baseline(&mut self) {
        let l = self.data.len();
        if l < 2 {
            return;
        }
        self.sort_by_x();
        let minmax = MinMax::new(self.data[0].y, self.data[l - 1].y);
        self.flatten_baseline_inner(&minmax);
    }

    /// A basic function to fill a point range with a straight line.
    /// NB, this function assumes a spectrum is x-sorted.
    pub fn fill_linear_gradient_within(&mut self, range: MinMax) {
        // early return clause. We do not work for silly spectra.
        if self.data.len() < 3 {
            return;
        }

        // If spectrum if long enough, find the minimum index for the gradient.
        let mut mindex = None;
        'loop1: for i in 1..self.data.len() {
            if self.data[i - 1].x == range.min {
                mindex = Some(i - 1);
                break 'loop1;
            } else if (self.data[i - 1].x < range.min) && (self.data[i].x > range.min) {
                mindex = Some(i);
                break 'loop1;
            }
        }

        let mindex = if let Some(index) = mindex {
            index
        } else {
            return;
        };

        // Next, find the maximum. NB, this is the index of the maximum
        // point used. When range is used, +1 should be used.
        let mut maxdex = None;
        'loop2: for i in 1..(self.data.len()) {
            if self.data[i - 1].x == range.max {
                maxdex = Some(i - 1);
                break 'loop2;
            } else if (self.data[i - 1].x < range.max) && (self.data[i].x > range.max) {
                maxdex = Some(i - 1);
                break 'loop2;
            }
        }

        let maxdex = if let Some(index) = maxdex {
            index
        } else {
            return;
        };

        // Get the step sizes for the gradient.
        let y_min = self.data[mindex].y;
        let y_max = self.data[maxdex].y;

        let dy = (y_max - y_min) / (maxdex + 1 - mindex) as f64;

        for i in mindex..(maxdex + 1) {
            // NB, with f64, this is more precise than multiple addition.
            let ans = y_min + dy * (i - mindex) as f64;
            //println!("point[i].y = {}", ans);
            self.data[i].y = ans;
        }
    }

    /// A function that finds peaks and then dedups them.
    /// Currently the width of the peak is not correct.
    pub fn find_peaks(&mut self, min_width: f64, max_width: f64, _min_height: f64, min_p: f64) {
        let ld = self.data.len();
        let x_range = self.data[ld - 1].x - self.data[0].x;
        let x_res = x_range / (ld - 1) as f64;

        let mut min_w = (min_width / x_res) as usize;
        if min_w < 1 {
            min_w = 1
        };
        let max_w = (max_width / x_res) as usize;

        // Determine the smallest acceptable peak.
        let min_y = self.data.iter().fold(
            std::f64::MAX,
            |acc, datum| if datum.y < acc { datum.y } else { acc },
        );

        let max_y = self.data.iter().fold(
            std::f64::MIN,
            |acc, datum| if datum.y > acc { datum.y } else { acc },
        );

        let min_p = (max_y - min_y) * min_p;

        // Iterate through each possible width value.
        for w in min_w..max_w {
            // iterate through the data, starting at points within a
            // "buffer" zone, to avoid going out of bounds.
            for i in (w..(ld - w)).rev() {
                let y = self.data[i].y;
                if self.data[(i - w)..(i + w + 1)].iter().fold(true, |acc, d| {
                    if y >= d.y {
                        acc
                    } else {
                        false
                    }
                }) {
                    let prominence = self.data[(i - w)..(i + w + 1)].iter().fold(0.0, |acc, d| {
                        let dy = y - d.y;
                        if dy > acc {
                            dy
                        } else {
                            acc
                        }
                    });

                    // If a peak is too small, we skip to the next point.
                    if prominence < min_p {
                        continue;
                    }
                    let peak = Peak::new(self.data[i].clone(), x_res * w as f64, prominence);

                    // Push only peaks that have a unique point.
                    // If a previously found peak is wider than previously found,
                    // push replace the width.
                    let mut dont_push = false;
                    for j in (0..self.peaks.len()).rev() {
                        if (self.peaks[j].max == peak.max) && (self.peaks[j].width < peak.width) {
                            dont_push = true;
                            self.peaks[j].width = peak.width;
                        } else if self.peaks[j].max == peak.max {
                            dont_push = true;
                        }
                    }
                    if !dont_push {
                        self.peaks.push(peak);
                    }
                }
            }
        }
        // The loop that makes peaks finishes here.
        // So now we can calculate half heights.
        for p in self.peaks.iter_mut() {
            let mut range = [None, None];
            let mut data = self.data.iter().peekable();
            let half_height = p.max.y / 2.0;

            // get the ranges.
            'loopa: while let Some(dat) = data.next() {
                if let Some(dat_b) = data.peek() {
                    if range[0].is_none() && (dat.y < half_height) && (dat_b.y >= half_height) {
                        let x =
                            dat.x + (dat_b.x - dat.x) * (half_height - dat.y) / (dat_b.y - dat.y);
                        range[0] = Some(x);
                    } else if range[0].is_some()
                        && (dat.y >= half_height)
                        && (dat_b.y < half_height)
                    {
                        let x = dat_b.x
                            - (dat_b.x - dat.x) * (half_height - dat_b.y) / (dat.y - dat_b.y);
                        range[1] = Some(x);
                        break 'loopa;
                    }
                } else {
                    break 'loopa;
                }
            }

            // Set the half-height width from this range.
            if range[0].is_some() && range[1].is_some() {
                p.set_half_h_width(range[1].unwrap() - range[0].unwrap());
            }
        }
    }

    /// A function to normalise the prominence of the peaks, setting the highest at 1.0.
    pub fn normalise_peak_prominence(&mut self) {
        // Get the maximum prominence of all peaks in the speactrum.
        let max_prominence = self.peaks.iter().fold(std::f64::MIN, |acc, p| {
            if p.prominence > acc {
                p.prominence
            } else {
                acc
            }
        });

        for peak in self.peaks.iter_mut() {
            peak.prominence /= max_prominence;
        }
    }

    /// Subtracts two spectra and creates a third.
    /// The third spectrum is shortened to the minimum range.
    /// NB, this is a heavy function which copies both spectra before
    /// working with them, in order to simplify things.
    pub fn subtract_spectrum(&self, other: &Spectrum) -> Option<Spectrum> {
        // If our spectra are two short to work with return a `None`.
        let a_len = self.data.len();
        let b_len = other.data.len();
        if (a_len < 3) || (b_len < 3) {
            return None;
        }

        // Clone and sort the data.
        let mut a = self.clone();
        a.sort_by_x();
        let mut b = other.clone();
        b.sort_by_x();

        let mut a_iter_mut = a.data.iter_mut().peekable();

        // Skip the points on the initial vector that are too low.
        'loop1: while let Some(dat) = a_iter_mut.peek() {
            if (b.data[0].x == dat.x) || ((b.data[0].x < dat.x) && (b.data[1].x > dat.x)) {
                break 'loop1;
            } else {
                a_iter_mut.next();
            }
        }

        // If there is no overlap between the two spectra then there
        // is nothing to subtract!
        if a_iter_mut.peek().is_none() {
            return None;
        }

        // Prepare struct for the new spectrum.
        let min_l = if a_len > b_len { b_len } else { a_len };
        let mut data = Vec::with_capacity(min_l);
        let peaks = Vec::new();

        // Work with the two spectra.
        let mut b_iter = b.data.iter().peekable();
        'loop2: while let Some(dat) = a_iter_mut.next() {
            if let Some(b_dat1) = b_iter.next() {
                // If by some miracle the x points match do a simple sub.
                //println!("dat={:?}\ndatb={:?}",dat,b_dat1);
                if dat.x == b_dat1.x {
                    let point = XYDatum::new(dat.x, dat.y - b_dat1.y);
                    data.push(point);
                } else if let Some(b_dat2) = b_iter.peek() {
                    // If this point is still not at the x-coordinate of
                    // `dat`, we need to keep going.
                    if b_dat2.x <= dat.x {
                        continue 'loop2;
                    // otherwise we interpolate between points.
                    } else {
                        let step = (dat.x - b_dat1.x) / (b_dat2.x - b_dat1.x);
                        let y_virtual = b_dat1.y + step * (b_dat2.y - b_dat1.y);
                        let point = XYDatum::new(dat.x, dat.y - y_virtual);
                        data.push(point);
                    }
                } else {
                    // NB: If we run out of points in spectrum_b we finish.
                    break 'loop2;
                }
            } else {
                // NB: If we run out of points in spectrum_b we finish.
                break 'loop2;
            }
        }

        Some(Spectrum { data, peaks })
    }

    /// Function to create a masked spectrum from your own initial spec,
    /// and a second containing peaks.
    pub fn create_masked_from(&self, other: &Spectrum) -> Spectrum {
        // If the masking spectrum has no peaks, then we do nothing.
        if other.peaks.is_empty() {
            return Spectrum::new();
        }
        //pub max: XYDatum,
        //pub width: f64,
        //pub prominence: f64,

        // Initially the output should contain the point data from self
        // and no peak data.
        let mut output = Spectrum::new();
        output.data = self.data.clone();

        // Iterate through each peak and mask it.
        for peak in other.peaks.iter() {
            let min = peak.max.x - peak.width / 2.0;
            let max = peak.max.x + peak.width / 2.0;
            let range = MinMax::new(min, max);
            output.fill_linear_gradient_within(range);
        }

        output
    }

    /// A function that tells you whether peaks are similar or not.
    /// This will only work when comparing samples with a single asbestos.
    pub fn compare_peaks_pure(
        &self,
        other: &Spectrum,
        no_peaks: usize,
        spectral_noise: f64,
    ) -> bool {
        // Filter out peaks that are too small to be of significance.
        // Then order by size.
        // println!("Own peaks:{:?}\nother peaks:{:?}",self.peaks,other.peaks);
        let mut own_sig_peaks = self
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .take(no_peaks)
            .collect::<Vec<_>>();
        own_sig_peaks.sort_by(|a, b| {
            a.prominence
                .partial_cmp(&b.prominence)
                .expect("Failed when sorting because peak prominence contained NaN.")
        });

        // Same here.
        let mut other_sig_peaks = other
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .take(no_peaks)
            .collect::<Vec<_>>();
        other_sig_peaks.sort_by(|a, b| {
            a.prominence
                .partial_cmp(&b.prominence)
                .expect("Failed when sorting because peak prominence contained NaN.")
        });

        // This is necessary.
        if own_sig_peaks.is_empty() || other_sig_peaks.is_empty() {
            return false;
        };

        // This is also necessary or the zip does not do its work.
        if own_sig_peaks.len() > other_sig_peaks.len() {
            return false;
        }

        // If the position of the peaks in order of size is the same, we say they match.
        for (own, other) in own_sig_peaks.iter().zip(other_sig_peaks.iter()) {
            //println!("Own peak: {:?}\nOther Peak: {:?}",own,other);
            let o_min_bound = other.max.x - X_TOLERANCE;
            let o_max_bound = other.max.x + X_TOLERANCE;
            if (o_max_bound < own.max.x) || (own.max.x < o_min_bound) {
                return false;
            }
        }

        true
    }

    /// A function that tells you whether peaks are similar or not.
    /// Does not order the peaks of the sample, so can analyse mixtures somewhat.
    pub fn compare_peaks_leniant(
        &self,
        other: &Spectrum,
        no_peaks: usize,
        tolerance: f64,
        spectral_noise: f64,
    ) -> bool {
        // Filter out peaks that are too small to be of significance.
        // Then order by size.
        // println!("Own peaks:{:?}\nother peaks:{:?}",self.peaks,other.peaks);
        let mut own_sig_peaks = self
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .collect::<Vec<_>>();
        own_sig_peaks.sort_by(|a, b| {
            a.prominence
                .partial_cmp(&b.prominence)
                .expect("Failed when sorting because peak prominence contained NaN.")
        });

        // We need to actually order peaks, biggest to smallest and
        // take only the most prominent.
        let own_sig_peaks = own_sig_peaks
            .iter()
            .rev()
            .take(no_peaks)
            .collect::<Vec<_>>();

        // Same here.
        let other_sig_peaks = other
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .collect::<Vec<_>>();

        // This is necessary.
        if own_sig_peaks.is_empty() || other_sig_peaks.is_empty() {
            return false;
        };

        // This is also necessary or the zip does not do its work.
        if own_sig_peaks.len() > other_sig_peaks.len() {
            return false;
        }

        // If the position of the peaks in order of size is the same, we say they match.
        for own in own_sig_peaks.iter() {
            //println!("Own peak: {:?}", own);
            let mut has = false;
            for other in other_sig_peaks.iter() {
                //println!("Own peak: {:?}\nOther Peak: {:?}",own,other);
                let o_min_bound = other.max.x - tolerance;
                let o_max_bound = other.max.x + tolerance;
                if (o_max_bound > own.max.x) && (own.max.x > o_min_bound) {
                    has = true;
                    break;
                }
            }
            if !has {
                return false;
            }
        }

        true
    }

    /// A function that tells you whether peaks are similar or not.
    /// Does not order sample peaks. Does measure x-distance between peaks.
    /// This is important for distinguishing crocidolite and amosite.
    pub fn compare_peaks_with_range(
        &self,
        other: &Spectrum,
        no_peaks: usize,
        spectral_noise: f64,
        tolerance: f64,
    ) -> bool {
        // Filter out peaks that are too small to be of significance.
        // Then order by size.
        // println!("Own peaks:{:?}\nother peaks:{:?}",self.peaks,other.peaks);
        let mut own_sig_peaks = self
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .collect::<Vec<_>>();
        own_sig_peaks.sort_by(|a, b| {
            a.prominence
                .partial_cmp(&b.prominence)
                .expect("Failed when sorting because peak prominence contained NaN.")
        });

        // We need to actually order peaks, biggest to smallest and
        // take only the most prominent.
        let own_sig_peaks = own_sig_peaks
            .iter()
            .rev()
            .take(no_peaks)
            .collect::<Vec<_>>();

        // Same here.
        let other_sig_peaks = other
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .collect::<Vec<_>>();

        // This is necessary.
        if own_sig_peaks.is_empty() || other_sig_peaks.is_empty() {
            return false;
        };

        // This is also necessary or the zip does not do its work.
        if own_sig_peaks.len() > other_sig_peaks.len() {
            return false;
        }

        // A variable that tells the function to examine distances between peaks
        // if there are at least two peaks.
        let do_range = own_sig_peaks.len() > 1;
        let mut other_x_points = [0.0; 2];

        // If the position of the peaks in order of size is the same, we say they match.
        for (i, own) in own_sig_peaks.iter().enumerate() {
            //println!("Own peak: {:?}", own);
            let mut has = false;
            for other in other_sig_peaks.iter() {
                //println!("Own peak: {:?}\nOther Peak: {:?}",own,other);
                let o_min_bound = other.max.x - tolerance;
                let o_max_bound = other.max.x + tolerance;
                if (o_max_bound > own.max.x) && (own.max.x > o_min_bound) {
                    if i < 2 {
                        other_x_points[i] = other.max.x;
                    }
                    has = true;
                    break;
                }
            }
            if !has {
                return false;
            }
        }

        // Test whether the range between the first two peaks is same
        // if applicable.
        if do_range {
            let own_rng = own_sig_peaks[0].max.x - own_sig_peaks[1].max.x;
            let oth_rng_min = other_x_points[0] - other_x_points[1] - tolerance;
            let oth_rng_max = other_x_points[0] - other_x_points[1] + tolerance;
            if (own_rng < oth_rng_min) || (oth_rng_max < own_rng) {
                return false;
            }
        }

        true
    }

    pub fn compare_peaks_ac_special(
        &self,
        other: &Spectrum,
        no_peaks: usize,
        tolerance: f64,
        spectral_noise: f64,
        min_mix_half_height_width: f64,
    ) -> (bool, bool) {
        // Filter out peaks that are too small to be of significance.
        // Then order by size.
        // println!("Own peaks:{:?}\nother peaks:{:?}",self.peaks,other.peaks);
        let mut own_sig_peaks = self
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .collect::<Vec<_>>();
        own_sig_peaks.sort_by(|a, b| {
            a.prominence
                .partial_cmp(&b.prominence)
                .expect("Failed when sorting because peak prominence contained NaN.")
        });

        // We need to actually order peaks, biggest to smallest and
        // take only the most prominent.
        let own_sig_peaks = own_sig_peaks
            .iter()
            .rev()
            .take(no_peaks)
            .collect::<Vec<_>>();

        // Same here.
        let other_sig_peaks = other
            .peaks
            .iter()
            .filter(|p| p.max.y > spectral_noise)
            .collect::<Vec<_>>();

        // This is necessary.
        if own_sig_peaks.is_empty() || other_sig_peaks.is_empty() {
            return (false, false);
        };

        // This is also necessary or the zip does not do its work.
        if own_sig_peaks.len() > other_sig_peaks.len() {
            return (false, false);
        }

        // A variable that tells the function to examine distances between peaks
        // if there are at least two peaks.
        let do_range = own_sig_peaks.len() > 1;
        let mut other_x_points = [0.0; 2];

        // If the position of the peaks in order of size is the same, we say they match.
        let mut other_peak_list = Vec::new();
        for (i, own) in own_sig_peaks.iter().enumerate() {
            //println!("Own peak: {:?}", own);
            let mut has = false;
            for (j, other) in other_sig_peaks.iter().enumerate() {
                //println!("Own peak: {:?}\nOther Peak: {:?}",own,other);
                let o_min_bound = other.max.x - tolerance;
                let o_max_bound = other.max.x + tolerance;
                if (o_max_bound > own.max.x) && (own.max.x > o_min_bound) {
                    other_peak_list.push(j);
                    if i < 2 {
                        other_x_points[i] = other.max.x;
                    }
                    has = true;
                    break;
                }
            }
            if !has {
                return (false, false);
            }
        }

        // We are assuming at least one peak is here, because we checked above.
        // If the main peak is too fat, we have both amosite and crocidolite.
        println!(
            "half_width is: {}",
            other_sig_peaks[other_peak_list[0]].half_h_width
        );
        if other_sig_peaks[other_peak_list[0]].half_h_width > min_mix_half_height_width {
            (true, true)
        // Test whether the range between the first two peaks is same
        // if applicable. This is an extra exclusion criteria.
        // If all else fails, we have either amosite or crocidolite.
        } else {
            if do_range {
                let own_rng = own_sig_peaks[0].max.x - own_sig_peaks[1].max.x;
                let oth_rng_min = other_x_points[0] - other_x_points[1] - tolerance;
                let oth_rng_max = other_x_points[0] - other_x_points[1] + tolerance;
                if (own_rng < oth_rng_min) || (oth_rng_max < own_rng) {
                    return (false, false);
                }
            }
            (true, false)
        }
    }
}
