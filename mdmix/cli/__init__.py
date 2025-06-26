from cli.base import app, DisplayUtils
from cli.project_plugin import ProjectPlugin
from cli.solvent_plugin import SolventPlugin

# Register plugins
project_plugin = ProjectPlugin()
project_plugin.attach_to_main(app)

solvent_plugin = SolventPlugin()
solvent_plugin.attach_to_main(app)

__all__ = ['app', 'DisplayUtils', 'project_plugin', 'solvent_plugin']