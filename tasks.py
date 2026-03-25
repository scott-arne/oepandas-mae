import sys
# noinspection PyPackageRequirements
from invoke.tasks import task
from pathlib import Path

ROOT = Path(__file__).parent.absolute()


@task
def test(c):
    c.run(f'cd {ROOT} && {sys.executable} -m unittest')

@task
def build(c):
    c.run(f'rm -rf {ROOT / "dist"}')
    c.run(f'cd {ROOT} && {sys.executable} -m build')

@task
def upload(c):
    c.run(f'scp {ROOT}/dist/* hopi36.pri.bms.com:/web/msdpypi/packages/ && echo "Successfully uploaded..."')

@task
def deploy(c):
    build(c)
    upload(c)

