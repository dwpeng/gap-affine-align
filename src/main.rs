use gap_affine_align::align;

fn main() {
    let text = "ATCGGATCTACTATCATCTACTA";
    let patt = "GATCTACTATCAT";
    let align_option = align::GapAffineAlignmentOptions {
        match_score: 2,
        mismatch_score: -4,
        gap_open: -2,
        gap_extend: -1,
    };
    let mut aligner = align::GapAffineAlignment::new(text.len(), patt.len(), align_option);
    let align_ret = aligner.align(text, patt);
    println!("{:?}", align_ret);
    align_ret.pretty_print(text, patt);
}
