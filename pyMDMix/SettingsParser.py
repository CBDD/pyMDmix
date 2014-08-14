##
## pyMDMix --- http://mdmix.sourceforge.net
## Software for preparation, analysis and quality control
## of solvent mixtures molecular dynamics
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
## Please cite your use of pyMDMix in published work:
##
##    TOBEPUBLISHED
##

"""
SettingsParser
==============

Read and parse ini style settings files. These type of files are used for project configuration,
replica configuration, solvent creation and settings parsing.

This module has been partly written from the original Biskit.SettingsParser and
Biskit.SettingsManager modules. Thanks to Raik.
"""

__author__="dalvarez"
__date__ ="$16-ene-2014 19:28:13$"


import logging
import os
import os.path as osp
import ConfigParser
import tools as T
import string

class SettingsError( Exception ):
    pass

class InvalidType( SettingsError ):
    pass

class InvalidValue( SettingsError ):
    pass

class InvalidFile( SettingsError ):
    pass

class SettingsWarning( SettingsError ):
    pass

class InvalidPath( SettingsWarning ):
    pass

class InvalidBinary( SettingsWarning ):
    pass

class WriteCfgError(SettingsError):
    pass


class CaseSensitiveConfigParser( ConfigParser.SafeConfigParser ):
    """
    Change ConfigParser so that it doesn't convert option names to lower case.
    """
    def optionxform(self, optionstr):
        return optionstr


class Setting(object):
    """
    Simple container for a single parameter
    """
    ## types of settings // section name in cfg file
    GENERAL = 'GENERAL'
    BIN = 'BINARIES'
    PATH = 'PATHS'

    def __init__( self, name=None, value=None, vtype=str, comment=None,
                  error=None, section=GENERAL ):
        self.name = name
        self.value = value
        self.vtype = vtype
        self.comment = comment
        self.error = error
        self.section = section

    def typeCast( self, vtype ):
        """
        Recast value to a new type. Value None remains unchanged.
        @param vtype: new type for value
        @type  vtype: type
        @raise InvalidValue: if current value is incompatible with vtype
        """
        try:
            if not self.value is None:
                if vtype == list and isinstance(self.value, str):
                    self.value = map(string.strip, self.value.split(','))
                else:
                    self.value = vtype( self.value )
            self.vtype = vtype
        except ValueError, e:
            raise InvalidValue, '%s: cannot convert "%s" to %r.' %\
              (self.name,self.value,vtype)

    def __repr__( self, tab='' ):
        error = ''
        if self.error: error = '(!)'

        comment = self.comment or ''
        if comment: comment = ' # '+comment

        return '%s%s = %s%s(%s)%s%s' %\
               (error, self.name, tab, self.vtype.__name__, str(self.value),\
               tab, comment)

    def __str__( self ):
        return self.__repr__( tab='\t' )

    def __cmp__( self, other ):
        """
        Compare Setting instances by their name.
        @return: 1,0,-1
        @rtype: int
        """
        if isinstance( other, self.__class__ ):
            return cmp( self.name, other.name )

        return cmp( self, other )


    def formatted( self ):
        """
        @return: parameter formatted for setting file
        @rtype: str
        """
        comment = ''
        error = ''
        name = self.name
        value = self.value or ''

        if self.vtype != str:
            name = self.vtype.__name__ + '-' + name
            if self.vtype == list: value = ','.join(value)

        if self.comment:
            comment = '\t## ' + self.comment

        if self.error:
            error = '\t#! ' + self.error + ' !'


        return '%s = %s%s%s' % (name, str(value), comment, error)


class SettingsParser(object):
    """
    A config file parser on steroids -- performs the following tasks:

      1. read a ini-style settings file
      2. type-cast options (e.g. of the form int-some_name into int(some_name))
      3. validate that all entries of section [PATHS] point to existing paths
      4. absolutize all valid paths
      5. validate that all entries of section [BINARIES] point to binaries
    """

    def __init__(self, ini):
        self.log = logging.getLogger("SettingsParser")
        self.f_ini = T.absfile( ini )
        self.result = {}

    def __type( self, option, default=str ):
        """
        Extract type from option name.

        @param option: name of parameter
        @type  option: str
        @param default: default type [str]
        @type  default: type

        @return: type, stripped option name (e.g. 'int_var1' -> int, 'var1')
        @rtype: type, str

        @raise TypeError: if type cannot be interpreted
        """
        t = default
        o = option

        if option.count('-') > 0:

            try:

                splt = option.split('-')

                s = splt[0]
                o = ''.join( splt[1:] )

                t = eval( s )

                if not type(t) is type:
                    raise TypeError, '%s is not a valid type' % s

            except Exception, e:
                raise TypeError, 'Cannot extract type from %s: %r'\
                      % option, e

        return t, o


    def __process( self, option, value, section=Setting.GENERAL ):
        """
        @param option: option name
        @type  option: str

        @param value: option value
        @type  value: str

        @param section: which section are we working on
        @type  section: str

        @return: new setting
        @rtype: Setting

        @raise SettingsError: InvalidType or Value
        """
        r = Setting( section=section )

        try:

            x = value.split('#')             ## split off comments
            r.value = x[0].strip() or None   ## don't return empty strings

            if len(x) > 1:
                r.comment = ''.join( x[1:] )

            vtype, r.name = self.__type( option )
            r.typeCast( vtype )

            if section == Setting.PATH:
                r.value = T.validPath( r.value )

            if section == Setting.BIN:
                r.value = T.validBinary( r.value )

        except SettingsWarning, e:           ## catch and record warnings
            r.error = str(e)

        return r


    def __processSection( self, items, section=Setting.GENERAL, verbose=False ):
        """
        @param items: section comming from ConfigParser
        @type  items: [ ( str, str ) ]

        @param section: which config section are we working on?
        @type  section: str

        @return: validated path values
        @rtype : dict, {str: Setting}
        """
        r = {}

        for name, value in items:

            s = self.__process( name, value, section )

            r[ s.name ] = s

            if verbose and s.error:
                self.log.warning(s.error)

        return r


    def parse( self, keepsections=False ):
        """
        @return: dict of type-cast params contained in fini
        @rtype: dict, {str: Setting}

        @raise IOError: if the settings file does not exist
        @raise SettingsError: (InvalidFile, InvalidValue, InvalidType)
        """
        try:
            ## read from file
            c = CaseSensitiveConfigParser()

            if c.read( self.f_ini ) != [ self.f_ini ]:
                raise IOError, 'Settings file %s not found.' % self.f_ini

            for section in c.sections():
                res = self.__processSection( c.items(section), section)
                if keepsections: self.result.update({section:res})
                else: self.result.update(res)

        except ConfigParser.Error, e:
            raise InvalidFile, 'Error parsing settings file %s: ' %\
                  self.f_ini + str(e)

        return self.result

    def __repr__( self):
        r = super( SettingsParser, self).__repr__()
        err = len( [ s for s in self.result.values() if s.error ] )
        r += ' -- %i entries, (!) %i errors' % (len( self.result ), err)

        values = self.result.values()
        values.sort()

        for v in values:
            r+= T.clipStr( '\n - %s' % str(v), 75)

        return r


class SettingsManager(object):
    """
    SettingsManager merges the parameters from a default and a user
    configuration file into a python module where they are published as
    normal fields. The general flow is like this::

      default.cfg ---[SettingsParser]---\
                                         [SettingsManager]--->[settings]
                                               /
                user.cfg---[SettingsParser]---/

    See L{pyMDMix.SettingsParser}
    See L{pyMDMix.settings}

    The default configurations should be located in:

    * C{pyMDMix/data/defaults/settings.cfg}      --> L{pyMDMix.settings}

    The user configurations are expected in files of the same name in
    C{~/.mdmix/}.
    """

    USER_HEADER = """
##     This is a pyMDMix user configuration file. The parameters in
##     this file are overriding the default parameters given in
##     %(fdefault)s.
##     If missing, pyMDMix creates a new user configuration file with
##     those parameters for which the default value seems
##     invalid. The remaining parameters are commented out.

##     Parameters in this file will be accessible from within python as
##     fields of Biskit.settings. For example::
##
##       leaprc = some/path/to/leaprc  # some comment
##
##     will lead to a variable in pyMDMix.settings::
##
##     >>> import pyMDMix.setting as S
##     >>> S.leaprc
##     >>> 'some/path/to/leaprc'

##     ...If, and only if, leaprc also exists in the default settings
##     file.  Parameters that are not listed in the default settings file
##     are ignored.

##     The default type of parameters is str. A prefix to the name like
##     'int-', 'float-', 'bool-', 'list-', etc. will be interpreted as
##     type-casting. For example::
##
##       float-nice_value = 10  # some comment
##
##     will lead to a variable in Biskit.settings::
##
##     >>> S.nice_value
##     >>> 10.0

##     Pay special attention to 'list-' type, as this type will try to chop
##     a comma separated string into chunks and build a list of it.

"""

    def __init__( self, fdefault, fuser, createmissing=False, verbose=1 ):
        """
        @param fdefault: default configuration file
        @type  fdedault: str
        @param fuser: user configuration file
        @type  fuser: str
        @param createmissing: create user config file if missing
        @type  createmissing: bool
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0
        """
        self.log = logging.getLogger("SettingsManager")
        self.verbose = verbose
        self.fdefault = fdefault
        self.fuser = fuser
        self.createmissing = createmissing
        self.fusermissing = not os.path.exists( T.absfile(fuser) )

        self.settings = []  #: will hold extracted Setting's

    def __update( self, cfg_default, cfg_user ):
        """
        Override default settings by valid (or equally invalid) user settings.

        :Arguments:
            @param *cfg_default*: settings read in from default file
            @type  *cfg_default*: dict {'str':SettingsParser.Setting}
            @param *cfg_user*   : settings read in from user config file
            @type  *cfg_user*   : dict {'str':SettingsParser.Setting}

        :Returns:
            @return: configuration with valid user settings overriding default ones
            @rtype: dict {'str':SettingsParser.Setting}
        """
        r = {}

        for name, default in cfg_default.items():

            next = cfg_user.get( name, default )

            if next.error > default.error:

                if self.verbose: self.log.warning(\
                    'User setting %s is reset to default (%r),\n\treason: %s'\
                    % (name, default.value, next.error)\
                    + '\n\tPlease check %s!' % self.fuser )

                next = default

            r[name] = next

        return r


    def collectSettings( self ):
        """
        Parse and combine default and user-defined config files.
        """
        try:
            pdefault = SettingsParser( self.fdefault )
            cdefault = pdefault.parse()

            try:
                puser = SettingsParser( self.fuser )
                cuser = puser.parse()

            except IOError, e:
                if self.verbose: self.log.warning(
                    'Could not find file with user-defined settings in %s' \
                    % self.fuser)

                cuser = {}

            self.settings = self.__update( cdefault, cuser )

        except SettingsError, e:
            self.log.error( str(e) )
            raise SettingsError, str(e)


    def writeUserSettings( self, errorsonly=False ):
        """
        Create a settings file with all options that are invalid with their
        default value.
        """
        try:
            T.backup( self.fuser )  ## create backup if file already exists

            fpath = os.path.dirname(self.fuser)
            if not os.path.exists( fpath ):
                if self.verbose:
                    self.log.warning('Creating folder %s for pyMDMix settings.'\
                                       %fpath )
                os.mkdir( fpath )

            sections = [Setting.GENERAL]
            r = {}

            for section in sections:

                r[ section ] = [ s for s in self.settings.values() \
                                 if s.section == section]
                r[ section ].sort()

            f = open( self.fuser, 'w' )

            f.write( SettingsManager.USER_HEADER % self.__dict__ )

            for section in sections:

                f.write( '[%s]\n' % section )
                f.write('\n')

                for param in r[section]:

                    if (not errorsonly) or param.error:
                        f.write( param.formatted() + '\n')
                    else:
                        f.write( '## ' + param.formatted() + '\n')

                f.write('\n')

            f.close()

        except OSError, e:
            raise WriteCfgError, e

    def settings2dict( self ):
        """
        Create dictionary from settings.
        @return: dictionary of parameter names (keys) and values
        @rtype: dict {str : any}
        """
        return dict( [ (s.name, s.value) for s in self.settings.values() ] )


    def updateNamespace( self, ns, keepdefined=True):
        """
        1. Parse in default configuration and user configuration file
        2. Merge the two, preferring valid user settings
        3. Create missing user configuration file if createmissing=True
        4. Insert parameters into the given namespace

        @param ns: namespace of a module ( obtained with locals() )
        @type  ns: dict {str:any}

        @param keepdefined: Keep unchanged variables already present in namespace.
        @type  keepdefined: Bool
        """
        self.collectSettings()

        if self.fusermissing and self.createmissing:
            if self.verbose:
                self.log.warning('Creating new user configuration file %s.' \
                                   % self.fuser)
            self.writeUserSettings( errorsonly=True )

        d = self.settings2dict()

        if keepdefined: ns.update((k,v) for (k,v) in d.iteritems() if k not in ns.keys())
        else: ns.update(d)


import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""
    def test_SettingsManager(self):
        """SettingsManager test"""

        f_in = osp.join(T.dataRoot(), 'defaults', 'settings.cfg')
        self.f_out =  osp.join(T.tempDir(), 'settings.cfg')

        self.m = SettingsManager( f_in, self.f_out,
                                  createmissing=True,
                                  verbose=self.local )

        ns = locals()             ## fetch local namespace

        self.m.updateNamespace( ns ) ## parse and insert options into namespace

        if self.local:
            globals().update( locals() ) ## publish namespace for debugging

        r = self.m.settings2dict()['testparam']

        self.assertEqual( r, 42) ## from 'int-testparam = 42' in settings.cfg

    def test_SettingsParser(self):
        """SettingsParser test"""
        self.f_out = ''
        f_in = osp.join(T.templatesRoot(),'solvent_template.cfg')
        parser = SettingsParser(f_in)
        print parser.parse()
        parser.result = {}
        print parser.parse(keepsections=True)
        
    def cleanUp(self):
        if self.f_out: T.tryRemove( self.f_out, tree=1 )


if __name__ == '__main__':
    print "testing"
    BT.localTest()
