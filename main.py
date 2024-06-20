from dG_SEM_flux import logger
from dG_SEM_flux import dg_elastic_physical_fluxes, dg_elastic_nonlinear_friction

STAGE_NAME = "ELASTIC WAVE SIMULATION"

try:
    logger.info(f">>>>>> stage {STAGE_NAME} started <<<<<<") 
    # do something here to run simulations
    logger.info(f">>>>>> stage {STAGE_NAME} completed <<<<<<\n\nx==========x")
except Exception as e:
   logger.exception(e)
   raise e

STAGE_NAME = "EARTHQUAKE RUPTURE DYNNAMICS"

try:
    logger.info(f">>>>>> stage {STAGE_NAME} started <<<<<<") 
    # do something here to run simulations
    logger.info(f">>>>>> stage {STAGE_NAME} completed <<<<<<\n\nx==========x")
except Exception as e:
   logger.exception(e)
   raise e