mod backward_facing_step;
pub mod lid_driven_cavity;
mod pipe_flow;
pub mod plotter;
mod tools;

use backward_facing_step::BackwardFacingStep;
use pipe_flow::PipeFlow;
use plotter::Env;

use crate::lid_driven_cavity::LidDrivenCavity;
use crate::plotter::Plot;

pub trait Case {
    fn iterate(&mut self);
    fn has_converged(&self, epsilon: f64) -> bool;
    fn has_diverged(&self) -> bool;
    fn postprocessing(&self, plot: &mut Plot, iter: u32);
}

pub enum Cases {
    LidDrivenCavity(LidDrivenCavity),
    PipeFlow(PipeFlow),
    BackwardFacingStep(BackwardFacingStep),
}

impl Cases {
    pub fn new(
        case: &str,
        nx: usize,
        ny: usize,
        re: f64,
        relax_uv: f64,
        relax_p: f64,
        damping: f64,
    ) -> Self {
        match case {
            "lid_driven_cavity" => {
                Cases::LidDrivenCavity(LidDrivenCavity::new(nx, ny, re, relax_uv, relax_p, damping))
            }
            "pipe_flow" => Cases::PipeFlow(PipeFlow::new(nx, ny, re, relax_uv, relax_p, damping)),
            "backward_facing_step" => Cases::BackwardFacingStep(BackwardFacingStep::new(
                nx, ny, re, relax_uv, relax_p, damping,
            )),
            _ => unreachable!(),
        }
    }

    pub fn iterate(&mut self) {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.iterate(),
            Cases::PipeFlow(pipe_flow) => pipe_flow.iterate(),
            Cases::BackwardFacingStep(backward_facing_step) => backward_facing_step.iterate(),
        }
    }

    pub fn has_converged(&self, epsilon: f64) -> bool {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.has_converged(epsilon),
            Cases::PipeFlow(pipe_flow) => pipe_flow.has_converged(epsilon),
            Cases::BackwardFacingStep(backward_facing_step) => {
                backward_facing_step.has_converged(epsilon)
            }
        }
    }

    pub fn has_diverged(&self) -> bool {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.has_diverged(),
            Cases::PipeFlow(pipe_flow) => pipe_flow.has_diverged(),
            Cases::BackwardFacingStep(backward_facing_step) => backward_facing_step.has_diverged(),
        }
    }

    pub fn postprocessing(&self, iter: u32) {
        let env = Env::new();
        let mut plot = Plot::new(&env);

        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => {
                lid_driven_cavity.postprocessing(&mut plot, iter)
            }
            Cases::PipeFlow(pipe_flow) => pipe_flow.postprocessing(&mut plot, iter),
            Cases::BackwardFacingStep(backward_facing_step) => {
                backward_facing_step.postprocessing(&mut plot, iter)
            }
        }
    }
}
