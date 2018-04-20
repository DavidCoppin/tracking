'''
Created on April 20, 2018

@author: David Coppin
@institution: Department of Physics, University of Auckland

@description: Manages the parameters that will be used by the algorithm. Takes some from config.py.
              The rest (time, lat and lon) is given by the command line and is set by default 
              somewhere (here?)
'''

import argparse
import numpy as np
import os,sys,string

class Configuration(dict):

    def __init__(self, filename):
        """
        Try to open the file filename and read it
        @filename: the name of the config-file
        """
        try:
            f=open(filename)
            file = f.readlines()
            f.close()
        except IOError:
            print "Could not open %s: No such file or directory"
            return
        self.__get(file)
        for key,value in self.iteritems():
            try:
                if "$" in value:
                    var=value.split('/')[0]
                    try:
                        path=os.environ[var.replace('$','')]
                        self[key]=value.replace(var,path)
                    except KeyError:
                        raise KeyError('Environment variable %s not set'%var)
                elif '~' in value:
                    self[key] = value.replace('~',os.environ['HOME'])
            except TypeError:
                pass


    def __setattr__(self,k,v):
        """
        Set class attribute
        """
        if k in self.keys():
            self[k] = v
        elif not hasattr(self,k):
            self[k] = v
        else:
            raise AttributeError, "Cannot set '%s', class attribute already exists" % (k, )

    def __getattr__(self, k):
        """
        Get attribute
        """
        if k in self.keys():
            return self[k]
        raise AttributeError, "Attribute '%s', deos not exist, available attirbutes are: %s" \
                %(k,self.keys().__repr__().strip(']').strip('['))

    def __repr__(self):
        a=8
        b=8
        for v,k in self.iteritems():
            if len(str(v)) > a:
                a = len(str(v))
            if len(str(k)) > b:
                b = len(str(k))
        a = a + 1
        b = b + 1
        s="Varnames"+(a-8)*' '+"| Values\n"
        s=s+(a+b)*'_'+'\n'
        for v,k in self.iteritems():
            s=s+str(v)+(a-len(str(v)))*' '+'| '+str(k)+'\n'
        return s

    def __get(self,file):
        """
        Take a variable-name and extracts the according value form the config file
        @file : the name of the file where the variables are
        """
        #Define a blacklist of characters
        blacklist=[']','[','{','}','@','#','"',"'"]
        for i,line in enumerate(file):
            try:
                if not line.strip('\n')[0] in blacklist and '=' in line:
                    if len(line.strip('\n').split('=')[0].split(' ')[:-1])<=1:
                        var=line.strip('\n').strip('\t').strip()\
                                .split('=')
                        value=var[1].replace(' ','')
                        if (( b in blacklist )for b in value):
                            pos = len(value)
                            for b in value:
                                if b in blacklist:
                                    if value.index(b) <= pos:
                                        pos = value.index(b)
                            value=value[:pos]
                            try:
                                value = int(value)
                            except ValueError:
                                try :
                                    value = float(value)
                                except ValueError:
                                    if value.lower() == 'false':
                                        value = False
                                    elif value.lower() == 'true':
                                        value = True
                                    elif value.lower() == 'none':
                                        value = None
                        self[var[0].replace(' ','')]=value
            except IndexError:
                pass


