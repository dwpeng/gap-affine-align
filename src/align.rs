#[derive(Debug, Clone)]
pub struct GapAffineAlignmentOptions {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapAffineAlignmentDirection {
    Init,
    EO,
    FO,
    HX,
    HM,
}

pub struct GapAffineAlignment {
    pub options: GapAffineAlignmentOptions,
    v: Vec<i32>,
    e: Vec<i32>,
    f: Vec<i32>,
    b: Vec<Vec<GapAffineAlignmentDirection>>,
    result: GapAffineAlignmentResult,
}

#[derive(Debug)]
pub struct GapAffineAlignmentResult {
    pub mismatch_base: i32,
    pub match_base: i32,
    pub deletion: i32,
    pub insertion: i32,
    pub alignment: Vec<(usize, usize, GapAffineAlignmentDirection)>,
}

impl GapAffineAlignmentResult {
    fn default() -> Self {
        Self {
            mismatch_base: 0,
            match_base: 0,
            deletion: 0,
            insertion: 0,
            alignment: Vec::with_capacity(300),
        }
    }
}

impl GapAffineAlignmentResult {
    pub fn pretty_print(&self, text: &str, patt: &str) {
        let mut ref_seq = String::new();
        let mut query_seq = String::new();
        let mut aln_seq = String::new();
        for (i, j, aln) in self.alignment.iter() {
            match aln {
                GapAffineAlignmentDirection::HM => {
                    ref_seq.push(text.chars().nth(*i).unwrap());
                    query_seq.push(patt.chars().nth(*j).unwrap());
                    aln_seq.push('|');
                }
                GapAffineAlignmentDirection::HX => {
                    ref_seq.push(text.chars().nth(*i).unwrap());
                    query_seq.push(patt.chars().nth(*j).unwrap());
                    aln_seq.push('*');
                }
                GapAffineAlignmentDirection::FO => {
                    ref_seq.push(text.chars().nth(*i).unwrap());
                    query_seq.push('-');
                    aln_seq.push('-');
                }
                GapAffineAlignmentDirection::EO => {
                    ref_seq.push('-');
                    query_seq.push(patt.chars().nth(*j).unwrap());
                    aln_seq.push('-');
                }
                _ => {}
            }
        }
        println!("{}", ref_seq);
        println!("{}", aln_seq);
        println!("{}", query_seq);
    }
}

impl GapAffineAlignment {
    pub fn new(reflen: usize, qrylen: usize, options: GapAffineAlignmentOptions) -> Self {
        let v = vec![0; reflen + 1];
        let e = vec![0; reflen + 1];
        let f = vec![0; reflen + 1];
        let b = vec![vec![GapAffineAlignmentDirection::Init; qrylen + 1]; reflen + 1];
        Self {
            options,
            v,
            e,
            f,
            b,
            result: GapAffineAlignmentResult::default(),
        }
    }

    fn init_matrix(&mut self, _text: &str, patt: &str) {
        // set initial values
        for i in 1..patt.len() + 1 {
            self.v[i] = self.options.gap_open + i as i32 * self.options.gap_extend;
        }
        // set e
        self.e[0] = i32::MIN / 2;
        for i in 1..patt.len() + 1 {
            self.e[i] = self.options.gap_open + i as i32 * self.options.gap_extend;
        }
        // set f
        for i in 0..patt.len() + 1 {
            self.f[i] = i32::MIN / 2;
        }
    }

    fn backtrace(&mut self, text: &str, patt: &str) {
        let mut i = text.len();
        let mut j = patt.len();
        self.result.alignment.clear();
        self.result.deletion = 0;
        self.result.insertion = 0;
        self.result.mismatch_base = 0;
        self.result.match_base = 0;
        while i > 0 && j > 0 {
            self.result.alignment.push((i - 1, j - 1, self.b[i][j]));
            match self.b[i][j] {
                GapAffineAlignmentDirection::HX => {
                    i -= 1;
                    j -= 1;
                    self.result.mismatch_base += 1;
                }
                GapAffineAlignmentDirection::HM => {
                    i -= 1;
                    j -= 1;
                    self.result.match_base += 1;
                }
                GapAffineAlignmentDirection::FO => {
                    i -= 1;
                    self.result.deletion += 1;
                }
                GapAffineAlignmentDirection::EO => {
                    j -= 1;
                    self.result.insertion += 1;
                }
                GapAffineAlignmentDirection::Init => {
                    panic!("Invalid backtrace direction");
                }
            }
        }

        while j > 0 {
            self.result
                .alignment
                .push((0, j - 1, GapAffineAlignmentDirection::EO));
            j -= 1;
            self.result.insertion += 1;
        }

        while i > 0{
            self.result.alignment.push((i - 1, 0, GapAffineAlignmentDirection::FO));
            i -= 1;
            self.result.deletion += 1;
        }

        self.result.alignment.reverse();
    }

    pub fn align(&mut self, text: &str, patt: &str) -> &GapAffineAlignmentResult {
        self.init_matrix(text, patt);
        let h = match self.options.gap_open > 0 {
            true => self.options.gap_open,
            false => -self.options.gap_open,
        };
        let s = match self.options.gap_extend > 0 {
            true => self.options.gap_extend,
            false => -self.options.gap_extend,
        };
        let mut vl: i32;
        let mut e: i32;
        let mut f: i32;
        let mut x: i32;
        for (ii, base1) in text.chars().enumerate() {
            let i = ii + 1;
            self.v[0] = -(h + ((i - 1) as i32) * s);
            vl = self.v[0];
            self.f[0] = self.v[0];

            if i == 1 {
                self.v[0] = 0;
                vl = 0;
                self.f[0] = i32::MIN / 2;
            }
            for (jj, base2) in patt.chars().enumerate() {
                let j = jj + 1;
                e = std::cmp::max(self.e[j - 1], vl - h - s);
                self.e[j] = e;
                f = std::cmp::max(self.f[j] - s, self.v[j] - h - s);
                self.f[j] = f;
                x = match base1 == base2 {
                    true => self.options.match_score,
                    false => self.options.mismatch_score,
                };
                x += self.v[j - 1];
                self.v[j - 1] = vl;
                vl = if x > e {
                    if x > f {
                        self.b[i][j] = match base1 == base2 {
                            true => GapAffineAlignmentDirection::HM,
                            false => GapAffineAlignmentDirection::HX,
                        };
                        x
                    } else {
                        self.b[i][j] = GapAffineAlignmentDirection::FO;
                        f
                    }
                } else {
                    if e > f {
                        self.b[i][j] = GapAffineAlignmentDirection::EO;
                        e
                    } else {
                        self.b[i][j] = GapAffineAlignmentDirection::FO;
                        f
                    }
                };
            }
            self.v[patt.len()] = vl;
        }
        self.backtrace(text, patt);
        &self.result
    }
}
