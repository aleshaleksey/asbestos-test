// This is a structure corresponding to an x-y datum.
#[derive(Debug, Clone, PartialEq)]
pub struct XYDatum {
    pub x: f64,
    pub y: f64,
}

//A min-max structure.
#[derive(Debug, Clone)]
pub struct MinMax {
    pub min: f64,
    pub max: f64,
}

/// This structure contains x-ranges to be exluded from baseline calcs.
#[derive(Debug, Clone)]
pub struct XRanges {
    pub ranges: Vec<MinMax>,
}

/// This structure describes a peak on a spectrum.
#[derive(Debug, Clone)]
pub struct Peak {
    // The XY point corresponding to the peak maximum.
    pub max: XYDatum,
    // The "width" of the peak at its base.
    pub width: f64,
    // Width at the half-height of the peak.
    pub half_h_width: f64,
    // The relative height of the peak compared to the tallest in the spectrum.
    pub prominence: f64,
}

impl MinMax {
    /// A function to make a new min-max.
    pub fn new(min: f64, max: f64) -> MinMax {
        MinMax { min, max }
    }
}

impl Peak {
    /// Construct a new instance of Peak based on point and width.
    /// NB: `half_h_width` is not used by this function and is added later.
    pub fn new(max: XYDatum, width: f64, prominence: f64) -> Peak {
        let half_h_width = 0.0;
        Peak {
            max,
            width,
            half_h_width,
            prominence,
        }
    }

    /// This needless function inserts the `half_h_width`.
    pub fn set_half_h_width(&mut self, hhw: f64) {
        self.half_h_width = hhw;
    }
}

impl XYDatum {
    /// Create a blank datum which displays default values.
    pub fn blank() -> XYDatum {
        XYDatum {
            x: Default::default(),
            y: Default::default(),
        }
    }

    /// Set the x value of a datum to be a given value.
    pub fn set_x(&mut self, x: f64) {
        self.x = x;
    }

    /// Set the y value for a datum to be a given value.
    pub fn set_y(&mut self, y: f64) {
        self.y = y;
    }

    /// Create a datum with a set x and y value.
    pub fn new(x: f64, y: f64) -> XYDatum {
        XYDatum { x, y }
    }

    /// Get the xy points for a datum as f64.
    pub fn get_xy(&self) -> (f64, f64) {
        (self.x, self.y)
    }
}
