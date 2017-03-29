from md_engine import MdEngine
import hoomd
from hoomd import deprecated
from hoomd import md


def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False


class HOOMDEngine(MdEngine):

    def __init__(self, dt=1e-2):
        MdEngine.__init__(self, 'HOOMD')
        self.dt = dt
        if run_from_ipython():
            print('Initializing HOOMD in ipython')
            hoomd.context.initialize('--mode=cpu')
        else:
            hoomd.context.initialize()
        print('HOOMDEngine initialized.')

    def set_initial_structure(self, initial_structure):
        hoomd.init.read_snapshot(initial_structure)

    @staticmethod
    def read_initial_struct_from_file(file_name):
        if file_name.endswith('.hoomdxml'):
            system = hoomd.deprecated.init.read_xml(file_name)
        elif file_name.endswith('.gsd'):
            system = hoomd.init.read_gsd(file_name)
        return system

    def run_initial_structure(self, mix_time, output_dir):
        md.integrate.mode_standard(dt=self.dt)
        md.integrate.nve(group=hoomd.group.all())
        hoomd.run(mix_time)
        deprecated.dump.xml(group=hoomd.group.all(), filename=output_dir + "mix.hoomdxml", all=True)

    def run_md(self, run_time, output_dir):
        # Now we bond!
        # BOND=False
        # if BOND is True:
        #    # dpd.set_params(kT = bond_kT)
        #    bond_callback = hoomd.analyze.callback(callback=find_pair, period=bond_period)
        #    temp_log = hoomd.analyze.log(filename=None, quantities=["temperature"], period=bond_period)
        #    hoomd.run(bond_time)
        #    bond_callback.disable()
        #    temp_log.disable()
        #    deprecated.dump.xml(group=hoomd.group.all(), filename=cwd + run_dir + "bond.hoomdxml", all=True)
        # Now we run to eql
        # dpd.set_params(kT = eql_kT)

        # NO HOLD

        hoomd.run(run_time)
        deprecated.dump.xml(group=hoomd.group.all(), filename=output_dir + "final.hoomdxml", all=True)
