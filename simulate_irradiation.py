# Creates a sphere model of material with a fusion relevant point source
# simulates the particle transport and activation/depletion/transmutation

import math
import openmc
import openmc.deplete
from elements_to_simulate import elements
from pathlib import Path

import openmc_depletion_plotter
# this package provides convenient plotting methods for depletion simulations like this one
# more details here https://github.com/fusion-energy/openmc_depletion_plotter

# sets nuclear data directories
# TODO put in some defaults here based on ENDF/B 8.0
# openmc.config[]

# TODO also add some fusion relevant materials e.g. Eurofer

Path.mkdir(Path('figs/atoms'), parents=True, exist_ok=True)
Path.mkdir(Path('figs/activity'), parents=True, exist_ok=True)

# loops through all the elements ignoring neutrons
for element in elements:

    mats = openmc.Materials()

    # makes a simple material
    my_material = openmc.Material() 
    my_material.add_element(element, 1, percent_type='ao')
    my_material.set_density('g/cm3', 1)  #TODO find the correct density


    sphere_radius = 100
    volume_of_sphere = (4/3) * math.pi * math.pow(sphere_radius, 3)
    my_material.volume = volume_of_sphere  # a volume is needed so openmc can find the number of atoms in the cell/material
    my_material.depletable = True  # depletable = True is needed to tell openmc to update the material with each time step

    materials = openmc.Materials([my_material])
    materials.export_to_xml()


    # creates sphere model
    # surfaces
    sph1 = openmc.Sphere(r=sphere_radius, boundary_type='vacuum')

    # cells, makes a simple sphere cell
    shield_cell = openmc.Cell(region=-sph1)
    shield_cell.fill = my_material

    # sets the geometry to the universe that contains just the one cell
    geometry = openmc.Geometry([shield_cell])

    # creates a 14.06MeV neutron point source
    source = openmc.Source()
    source.space = openmc.stats.Point((0, 0, 0))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.06e6], [1]) # TODO allow different sources
    source.particles = 'neutron'

    # Instantiate a Settings object
    settings = openmc.Settings()
    settings.batches = 2  # TODO this to be user arg
    settings.inactive = 0
    settings.particles = 10000  # TODO this to be user arg
    settings.source = source
    settings.run_mode = 'fixed source'

    # TODO create dose tally on surface

    model = openmc.model.Model(geometry, materials, settings)

    operator = openmc.deplete.CoupledOperator(
        model=model,
        normalization_mode="source-rate",  # set for fixed source simulation, otherwise defaults to fission simulation
        dilute_initial=0,  # set to zero to avoid adding small amounts of isotopes, defaults to adding small amounts of fissionable isotopes
        reduce_chain=True,  # reduced to only the isotopes present in depletable materials and their possible progeny
        reduce_chain_level=5,
    )

    # We define timesteps together with the source rate to make it clearer
    # TODO find fusion relevant time scales
    timesteps_and_source_rates = [
        (24, 1e20),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
        (24, 0),
    ]

    # Uses list Python comprehension to get the timesteps and source_rates separately
    timesteps = [item[0] for item in timesteps_and_source_rates]
    source_rates = [item[1] for item in timesteps_and_source_rates]


    # PredictorIntegrator has been selected as the depletion operator for this example as it is a fast first order Integrator
    # OpenMC offers several time-integration algorithms https://docs.openmc.org/en/stable/pythonapi/deplete.html#primary-api\n",
    # CF4Integrator should normally be selected as it appears to be the most accurate https://dspace.mit.edu/handle/1721.1/113721\n",
    integrator = openmc.deplete.PredictorIntegrator(
        operator=operator,
        timesteps=timesteps,
        source_rates=source_rates
    )

    integrator.integrate()

    results = openmc.deplete.ResultsList.from_hdf5("depletion_results.h5")

    atoms_fig = results.plot_atoms_vs_time(excluded_material=my_material, plotting_backend='matplotlib')
    atoms_fig.savefig(f'figs/atoms/element_{element}.png')
    activity_fig = results.plot_activity_vs_time(excluded_material=my_material, plotting_backend='matplotlib')
    activity_fig.savefig(f'figs/activity/element_{element}.png')

    # TODO plots dose

    # TODO plots decay heat
