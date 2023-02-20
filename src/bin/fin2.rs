use bogay::Newton;

fn main() {
    let points = vec![(8.1, 16.94410), (8.3, 17.56492), (8.7, 18.82091)];
    let y = Newton::new(points).eval(8.4);
    println!("{y}");
}
