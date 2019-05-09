extern crate ngslib;
extern crate bio;
extern crate clap;
extern crate rust_htslib;
extern crate bio_types;
extern crate itertools;
extern crate rayon;

use bio::io::bed;
use bio::io::bed::*;
use bio::io::bed::Record;
use bio_types::strand::*;
use bio_types::annot::contig::*;
use bio_types::annot::loc::*;
use rust_htslib::bam;
use ngslib::ngslibrary::*;
use ngslib::locus::shift::Shift;
use rust_htslib::prelude::*;
use itertools::Itertools;
use rayon::prelude::*;



fn main() {
    use clap::{App, Arg, ArgGroup};

    let matches = App::new("cutsites")
                          .version("0.1.0-alpha.1")
                          .author("Matt Lawlor <matt.a.lawlor@gmail.com>")
                          .about("Get cut site profiles.")
                          .arg(Arg::with_name("IBAM")
                               .help("indexed bam file")
                               .required(true)
                               .index(1))
                          .arg(Arg::with_name("BED")
                               .help("regions in bed format")
                               .required(true)
                               .index(2))
                          .arg(Arg::with_name("FLANK")
                               .help("bp to add to either side of center of specified bed regions")
                               .short("f")
                                .long("flank")
                                .takes_value(true))
                          .arg(Arg::with_name("NAME")
                               .short("n")
                               .long("name")
                               .help("sample name")
                               .takes_value(true))
                          .arg(Arg::with_name("THREADS")
                          	   .help("threads to use")
                          	   .short("p")
                          	   .long("threads")
                          	   .takes_value(true))
                          .get_matches();

    let bed_file: &str = matches.value_of("BED").unwrap();
    let bam_file: &str = matches.value_of("IBAM").unwrap();
    let sample: &str = matches.value_of("NAME").unwrap_or(bam_file);
    let flank: usize = matches.value_of("FLANK").unwrap_or("50").parse().unwrap();

    let threads: usize = matches.value_of("THREADS").unwrap_or("1").parse().unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    eprintln!("Reading in bed regions...");
    let ctigs = bed_to_contigs(bed_file, flank);

    ctigs.par_chunks(1000)
         .for_each(|chunk| {
             let amap = bam_to_annotmap(bam_file, chunk);
             for c in chunk.iter() {
                 let cov = cov_across_contig(&c, &amap);
                 let r = format!("{}", c);
                 for (i,v) in cov.iter().enumerate() {
                     println!("{},{},{},{}", sample, r, i, v);
                 }
             }
         })

}

fn bed_to_contigs(b: &str,flank: usize) -> Vec<Contig<String,ReqStrand>> {
    let mut br = bed::Reader::from_file(b).unwrap();
    let res = br.records()
                .map(|a|a.unwrap())
                .map(|a| {
                    Contig::new(a.chrom().to_string(),
                        ((a.start() + a.end())/2) as isize - (flank as isize),
                        (flank * 2) as usize,
                        ReqStrand::Forward)});
                //.map(|a| pos_strandify(&a));
    res.collect()
}

fn tn5shift(c: Contig<String,ReqStrand>) -> Contig<String,ReqStrand> {
    match c.strand() {
        ReqStrand::Forward => c.shift(3).first_pos().contig(),
        ReqStrand::Reverse => c.shift(4).first_pos().contig(),
    }
}

fn bam_to_annotmap(bp: &str, cv: &[Contig<String,ReqStrand>]) -> NGSLibrary<bam::IndexedReader> {
    let mut bam = bam::IndexedReader::from_path(bp).unwrap();

    NGSLibrary::from_indexed(bam,
                            cv.to_vec(),
                            LibraryType::Unstranded,
                            Some(tn5shift))

}


fn cov_across_contig<T: Read>(c: &Contig<String,ReqStrand>, nl: &NGSLibrary<T>) -> Vec<usize>{
    nl.coverage_across(c)
}
