from epoxpy.bonding import FreudBonding
import hoomd
from hoomd import md
from hoomd import deprecated
import os
import numpy as np
from collections import Counter
from epoxpy.abc_type_epoxy_simulation import ABCTypeEpoxySimulation

class ResumeABCTypeEpoxySimulation(ABCTypeEpoxySimulation):
    """Simulation class for resuming an existing ABCTypeEpoxySimulation. This
    simulation uses the 'final.hoomdxml' to resume.
    """
    def set_initial_structure(self):
        self.initialize_system_from_file(self.resume_file_name)

    def initialize(self):
        print('Initializing resume for {}'.format(self.simulation_name))
        if not os.path.exists(self.output_dir):
            print('Creating simulation folder: {}'.format(self.output_dir))
            os.makedirs(self.output_dir)

        self.initialize_context()
        if self.ext_init_struct_path is None:
            self.set_initial_structure()
        else:
            print('Loading external initial structure : ', self.ext_init_struct_path)
            self.initialize_system_from_file(self.ext_init_struct_path)
            deprecated.dump.xml(group=hoomd.group.all(), filename=self.resume_file_name, all=True)

    def run(self):
        del self.system
        self.initialize_context()
        use_time_step_from_file = True
        self.initialize_system_from_file(self.resume_file_name, use_time_step_from_file=use_time_step_from_file)

        print('Running MD for {}'.format(self.simulation_name))
        self.setup_md_run()
        self.configure_outputs()
        if self.bond is True:
            log = hoomd.analyze.log(filename=None, quantities=["temperature"], period=self.bond_period)
            msd_groups = self.get_msd_groups()
            if self.legacy_bonding is True:
                self.bonding = LegacyBonding(system=self.system, groups=msd_groups, log=log,
                                                 activation_energy=self.activation_energy,
                                             sec_bond_weight=self.sec_bond_weight)
            else:
                self.bonding = FreudBonding(system=self.system, groups=msd_groups, log=log,
                                                activation_energy=self.activation_energy,
                                            sec_bond_weight=self.sec_bond_weight)
            bond_callback = hoomd.analyze.callback(callback=self.bonding, period=self.bond_period)

            if self.log_curing is True:
                curing_callback = hoomd.analyze.callback(callback=self.calculate_curing_percentage,
                                                         period=self.curing_log_period, phase=-1)
        hoomd.util.quiet_status()
        self.run_md()

        if self.bond is True:
            bond_callback.disable()
            log.disable()
            if self.log_curing is True:
                curing_callback.disable()

