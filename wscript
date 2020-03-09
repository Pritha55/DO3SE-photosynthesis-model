#!/usr/bin/env python3
import os
from waflib import Utils

APPNAME = 'fortran-environment'


if 'PREFIX' not in os.environ:
    os.environ['PREFIX'] = '.'

# Defines WAF options - https://waf.io/apidocs/Options.html
def options(ctx):
    # load compiler libraries (Adds commands to ctx)
    # https://waf.io/apidocs/tools/c_config.html
    ctx.load('compiler_c')
    # https://waf.io/apidocs/tools/fc_config.html
    ctx.load('compiler_fc')

def configure(ctx):
    print('Configuring')

    # load compilers and use them to check fortran
    ctx.load('compiler_c')
    ctx.load('compiler_fc')

    # Compiles a sample Fortran program to ensure that the settings are correct
    ctx.check_fortran()

def build(ctx):
    print('Building')
    print('Build Mocks')
    ctx(features='fc',
        source='src/mocked_functions.F90',
        target='modules',
        use='LAPACK C_API')

    print('Build B')

    ctx(features='fc',
        source=ctx.path.ant_glob('src/*_ml.F90'),
        target='modules',
        use='LAPACK C_API')

    ctx(features='fc',
        source=ctx.path.ant_glob('src/*_ml.f90'),
        target='modules',
        use='LAPACK C_API')

    print('Build C')

    ctx(features='fc fcprogram',
        source='src/main.f90',
        target='run_model',
        use='modules')