from dG_SEM_flux.constants import *
from dG_SEM_flux.utils.common import read_yaml, create_directories
from dG_SEM_flux.entity.config_entity import dGElasticConfig


class ConfigurationManager:
    def __init__(
        self,
        config_filepath = CONFIG_FILE_PATH,
        params_filepath = PARAMS_FILE_PATH):

        self.config = read_yaml(config_filepath)
        self.params = read_yaml(params_filepath)

        create_directories([self.config.artifacts_root])

    
    def get_dg_elastic_config(self) -> dGElasticConfig:
        config =  self.config.dGElastic
        params = self.params.dGElastic1D

        create_directories([config.root_dir])

        dg_elastic_config = dGElasticConfig(
            root_dir=config.root_dir,
            node=params.node,
            poly_deg=params.poly_deg,
            time_integrator=params.time_integrator
        )
        return dg_elastic_config