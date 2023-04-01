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
    fn iterate(&mut self) -> bool;
    fn postprocessing(&self, plot: &mut Plot, iter: usize);
}

pub enum Cases {
    LidDrivenCavity(LidDrivenCavity),
    PipeFlow(PipeFlow),
    BackwardFacingStep(BackwardFacingStep),
}

impl Cases {
    pub fn new(case: &str, nx: usize, ny: usize, re: f64) -> Self {
        match case {
            "lid_driven_cavity" => Cases::LidDrivenCavity(LidDrivenCavity::new(nx, ny, re)),
            "pipe_flow" => Cases::PipeFlow(PipeFlow::new(nx, ny, re)),
            "backward_facing_step" => {
                Cases::BackwardFacingStep(BackwardFacingStep::new(nx, ny, re))
            }
            _ => unreachable!(),
        }
    }

    pub fn postprocessing(&self, iter: usize) {
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

    pub fn iterate(&mut self) -> bool {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.iterate(),
            Cases::PipeFlow(pipe_flow) => pipe_flow.iterate(),
            Cases::BackwardFacingStep(backward_facing_step) => backward_facing_step.iterate(),
        }
    }
}
