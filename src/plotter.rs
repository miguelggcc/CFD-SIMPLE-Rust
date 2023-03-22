#![allow(dead_code)]
use std::path::Path;

use cpython::FromPyObject;
use cpython::NoArgs;
use cpython::ObjectProtocol; //for call method
use cpython::PyDict;
use cpython::PyList;
use cpython::PyObject;
use cpython::PyString;
use cpython::Python;
use cpython::PythonObject;
use cpython::{PyModule, PyTuple};

pub struct Env {
    gil: cpython::GILGuard,
}

impl Env {
    pub fn new() -> Env {
        Env {
            gil: Python::acquire_gil(),
        }
    }
}

pub struct Plot<'p> {
    py: Python<'p>,
    plt: PyModule,
    np: PyModule,
    artists: Vec<PyObject>,
}

impl<'p> Plot<'p> {
    pub fn new<'a>(env: &'a Env) -> Plot<'a> {
        let py = env.gil.python();
        let plt = PyModule::import(py, "matplotlib.pyplot").unwrap();
        let np = PyModule::import(py, "numpy").unwrap();

        let artists = Vec::default();
        Plot {
            py,
            plt,
            np,
            artists,
        }
    }

    pub fn ion(&self) {
        let _ = self.plt.call(self.py, "ion", NoArgs, None).unwrap();
    }

    pub fn ioff(&self) {
        let _ = self.plt.call(self.py, "ioff", NoArgs, None).unwrap();
    }

    pub fn show(&self) {
        let _ = self.plt.call(self.py, "show", NoArgs, None).unwrap();
    }

    pub fn imshow(&self, values: PyObject, cmap: &str, interpolation: &str) {
        let dict = PyDict::new(self.py);
        dict.set_item(self.py, "cmap", cmap).unwrap();
        dict.set_item(self.py, "interpolation", interpolation)
            .unwrap();
        let _ = self
            .plt
            .call(self.py, "imshow", (values,), Some(&dict))
            .unwrap();
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) {
        let dict = PyDict::new(self.py);
        dict.set_item(self.py, "bbox_inches", "tight").unwrap();
        let path = path.as_ref();
        let _ = self
            .plt
            .call(self.py, "savefig", (path.to_str().unwrap(),), Some(&dict))
            .unwrap();
    }

    pub fn scatter(&self, x: &[f32], y: &[f32]) {
        assert!(x.len() == y.len());
        let _ = self.plt.call(self.py, "scatter", (x, y), None).unwrap();
    }

    pub fn plot(&self, x: &[f32], y: &[f32]) {
        assert!(x.len() == y.len());
        let _ = self.plt.call(self.py, "plot", (x, y), None).unwrap();
    }

    pub fn semilogx(&self, x: &[f32], y: &[f32]) {
        assert!(x.len() == y.len());
        let _ = self.plt.call(self.py, "semilogx", (x, y), None).unwrap();
    }

    pub fn semilogy(&self, x: &[f32], y: &[f32]) {
        assert!(x.len() == y.len());
        let _ = self.plt.call(self.py, "semilogy", (x, y), None).unwrap();
    }

    pub fn pcolormesh(&mut self, x: &[f32], y: &[f32], z: &[f32], cmap: &str, title: &str) {
        let subplots = PyDict::new(self.py);
        subplots.set_item(self.py, "figsize", (6, 6)).unwrap();
        subplots.set_item(self.py, "dpi", 100).unwrap();
        let fig = self
            .plt
            .call(self.py, "figure", NoArgs, Some(&subplots))
            .unwrap();
        let ax = fig
            .call_method(self.py, "add_subplot", (111,), None)
            .unwrap();

        let dict = PyDict::new(self.py);
        dict.set_item(self.py, "cmap", cmap).unwrap();
        dict.set_item(self.py, "animated", "True").unwrap();

        dict.set_item(self.py, "vmin", 0.0).unwrap();
        dict.set_item(self.py, "vmax", 1.0).unwrap();

        let c = ax
            .call_method(
                self.py,
                "pcolorfast",
                (x, y, self.reshape(z, x.len(), y.len())),
                Some(&dict),
            )
            .unwrap();

        ax.call_method(self.py, "set_title", (title,), None)
            .unwrap();
        ax.call_method(self.py, "set_aspect", ("equal",), None)
            .unwrap();

        let dict2 = PyDict::new(self.py);
        dict2.set_item(self.py, "ax", ax).unwrap();
        dict2.set_item(self.py, "fraction", 0.046).unwrap();
        dict2.set_item(self.py, "pad", 0.04).unwrap();
        fig.call_method(self.py, "tight_layout", NoArgs, None)
            .unwrap();
        fig.call_method(self.py, "colorbar", (&c,), Some(&dict2))
            .unwrap();

        self.artists.push(c);
    }

    pub fn quiver(&mut self, x: &[f32], y: &[f32], u: &[f32], v: &[f32]) {
        let np_u = self.reshape(u, x.len(), y.len());
        let np_v = self.reshape(v, x.len(), y.len());

        let dictq = PyDict::new(self.py);
        dictq.set_item(self.py, "scale", 1 * x.len() / 30).unwrap();

        let q = self
            .plt
            .call(self.py, "quiver", (x, y, np_u, np_v), Some(&dictq))
            .unwrap();

        self.artists.push(q);
    }

    pub fn streamplot(&mut self, x: &[f32], y: &[f32], u: &[f32], v: &[f32]) {
        let np_u = self.reshape(u, x.len(), y.len());
        let np_v = self.reshape(v, x.len(), y.len());

        let dictq = PyDict::new(self.py);
        dictq.set_item(self.py, "density", 3.0).unwrap();

        let list = PyList::extract(
            self.py,
            &self.np.call(self.py, "meshgrid", (x, y), None).unwrap(),
        )
        .unwrap();

        let np_x = list.get_item(self.py, 0);
        let np_y = list.get_item(self.py, 1);

        self.plt
            .call(
                self.py,
                "streamplot",
                (np_x, np_y, np_u, np_v),
                Some(&dictq),
            )
            .unwrap();
    }

    pub fn contourf(&mut self, x: &[f32], y: &[f32], z: &[f32], cmap: &str, title: &str) {
        let subplots = PyDict::new(self.py);
        subplots.set_item(self.py, "figsize", (6, 6)).unwrap();
        subplots.set_item(self.py, "dpi", 100).unwrap();
        let fig = self
            .plt
            .call(self.py, "figure", NoArgs, Some(&subplots))
            .unwrap();
        let ax = fig
            .call_method(self.py, "add_subplot", (111,), None)
            .unwrap();

        let dict = PyDict::new(self.py);
        dict.set_item(self.py, "cmap", cmap).unwrap();

        let z_np = self.reshape(z, x.len(), y.len());

        let c = ax
            .call_method(self.py, "contourf", (x, y, z_np, 30), Some(&dict))
            .unwrap();

        ax.call_method(self.py, "set_title", (title,), None)
            .unwrap();
        ax.call_method(self.py, "set_aspect", ("equal",), None)
            .unwrap();

        let dict2 = PyDict::new(self.py);
        dict2.set_item(self.py, "ax", ax).unwrap();
        dict2.set_item(self.py, "fraction", 0.046).unwrap();
        dict2.set_item(self.py, "pad", 0.04).unwrap();
        fig.call_method(self.py, "tight_layout", NoArgs, None)
            .unwrap();
        fig.call_method(self.py, "colorbar", (&c,), Some(&dict2))
            .unwrap();
    }

    pub fn title(&self, title: &str) {
        let _ = self.plt.call(self.py, "title", (title,), None).unwrap();
    }

    pub fn xlabel(&self, label: &str) {
        let _ = self.plt.call(self.py, "xlabel", (label,), None).unwrap();
    }

    pub fn ylabel(&self, label: &str) {
        let _ = self.plt.call(self.py, "ylabel", (label,), None).unwrap();
    }

    pub fn grid(&self, grid: bool) {
        let _ = self.plt.call(self.py, "grid", (grid,), None).unwrap();
    }

    pub fn axis(&self, par: &str) {
        self.plt.call(self.py, "axis", (par,), None);
    }

    pub fn legend(&self, labels: &[&str]) {
        let gca = self.plt.call(self.py, "gca", NoArgs, None).unwrap();
        let labels: Vec<PyObject> = labels
            .iter()
            .map(|l| PyString::new(self.py, *l).into_object())
            .collect();
        let labels = PyList::new(self.py, &labels);
        gca.call_method(self.py, "legend", (labels,), None).unwrap();
    }

    pub fn draw(&self) {
        let _ = self
            .plt
            .call(self.py, "draw", PyTuple::empty(self.py), None)
            .unwrap();
    }

    pub fn clf(&self) {
        let _ = self
            .plt
            .call(self.py, "clf", PyTuple::empty(self.py), None)
            .unwrap();
    }

    pub fn reshape(&self, values: &[f32], a: usize, b: usize) -> PyObject {
        self.np
            .call(self.py, "reshape", (values, (a, b)), None)
            .unwrap()
    }
}

pub struct Animation<'p> {
    py: Python<'p>,
    plt: PyModule,
    plot: Plot<'p>,
    np: PyModule,
    writer: PyObject,
    artists: Vec<PyObject>,
}

impl<'p> Animation<'p> {
    pub fn new<'a>(env: &'p Env, plot: Plot<'p>, fps: u32) -> Animation<'p> {
        let py = env.gil.python();
        let plt = PyModule::import(py, "matplotlib.pyplot").unwrap();
        let np = PyModule::import(py, "numpy").unwrap();
        let animation = PyModule::import(py, "matplotlib.animation").unwrap();
        let dict = PyDict::new(py);
        dict.set_item(py, "fps", fps).unwrap();
        dict.set_item(py, "extra_args", vec!["-vcodec", "h264_nvenc"])
            .unwrap();
        let writer = animation
            .call(py, "FFMpegWriter", NoArgs, Some(&dict))
            .unwrap();
        let artists = Vec::default();
        Self {
            py,
            plot,
            plt,
            np,
            writer,
            artists,
        }
    }

    pub fn setup_animation(&self, filename: &str) {
        let fig = self.plt.call(self.py, "gcf", NoArgs, None).unwrap();
        self.writer
            .call_method(self.py, "setup", (fig, filename), None)
            .unwrap();
    }

    pub fn update_frame(&self, z: &[f32], _u: &[f32], _v: &[f32], x_len: usize, y_len: usize) {
        let np_z = self.reshape(z, x_len, y_len);
        self.artists[0]
            .call_method(self.py, "set_data", (&np_z,), None)
            .unwrap();

        self.writer
            .call_method(self.py, "grab_frame", NoArgs, None)
            .unwrap();
    }

    pub fn finish_animation(&self) {
        self.writer
            .call_method(self.py, "finish", NoArgs, None)
            .unwrap();
    }

    pub fn clf(&self) {
        let _ = self
            .plt
            .call(self.py, "clf", PyTuple::empty(self.py), None)
            .unwrap();
    }

    pub fn reshape(&self, values: &[f32], a: usize, b: usize) -> PyObject {
        self.np
            .call(self.py, "reshape", (values, (a, b)), None)
            .unwrap()
    }
}
