# dependency package for singulary-opac

import os
from spack import *

class SingularityOpacDeps(BundlePackage):
    homepage    = "https://github.com/lanl/singularity-opac"
    git         = "git@github.com:lanl/singularity-opac.git"

    version("main", branch="main")

    variant("cuda", default=False, description="Enable cuda support")

    depends_on("cmake", type="build")
    depends_on("hdf5~mpi+cxx+hl", type=("build", "run"))
    depends_on("kokkos@3:~shared+cuda+wrapper+serial+cuda_lambda+cuda_relocatable_device_code cuda_arch=70", when="+cuda", type=("build", "run"))
    depends_on("catch2", type=("build", "run"))

    phases=["install"]

    def setup_run_environment(self, env):
        env.set('HDF5_ROOT', self.spec['hdf5'].prefix)
        # this is still WIP
        env.set('SINGULARITY_EOS_LOADER', os.path.join(self.spec.prefix, f"load_env-{self.spec.full_hash(length=4)}.sh"))

    def install(self, spec, prefix):
        mod_script = os.path.join(spec.prefix, f"load_env-{spec.full_hash(length=4)}.sh")

        with open(os.path.join(mod_script), "w") as f:
            f.write(f"# load env {spec.short_spec}")
            f.write("")
            for dep in spec.dependencies(deptype="build"):
                f.write(dep.format())


