class apt-update {
    exec { "/usr/bin/apt-get -y update": }
}

class python {
    package {
        "build-essential": ensure => latest;
        "python": ensure => latest;
        "python-dev": ensure => latest;
        "python-setuptools": ensure => latest;
    }
    exec { "easy_install pip":
        path => "/usr/local/bin:/usr/bin:/bin",
        refreshonly => true,
        require => Package["python-setuptools"],
        subscribe => Package["python-setuptools"];
    }
}

class perl {
    package {
        "perl": ensure => latest;
    }
}

class python-libs {
    # ubuntu packages required to build the following python packages
    package {
        "libatlas-base-dev":
            ensure => "latest";
        "gfortran":
            ensure => "latest";
        "libpng12-dev":
            ensure => "latest";
        "python-gtk2-dev":
            ensure => "latest";
        "libfreetype6-dev":
            ensure => "latest";
        "libgraphviz-dev":
            ensure => "latest";
    }
    # python packages
    package {
        "numpy":
            ensure => "1.6.1",
            provider => "pip";
        "scipy":
            ensure => "0.10.1",
            provider => "pip",
            require => [
                Package["gfortran"],
                Package["libatlas-base-dev"],
                Package["numpy"]
            ];
        "sympy":
            ensure => "0.7.1",
            provider => "pip";
        "matplotlib":
            ensure => "1.1.0",
            provider => "pip",
            require => [
                Package["libpng12-dev"],
                Package["python-gtk2-dev"],
                Package["libfreetype6-dev"],
                Package["numpy"]
            ];
        "ipython":
            ensure => "0.12",
            provider => "pip";
        "pygraphviz":
            ensure => "1.1",
            provider => "pip",
            require => Package["libgraphviz-dev"];
    }
}

class bionetgen {
    package { "wget": ensure => present; }
    exec { "/usr/bin/wget http://dl.dropbox.com/u/19644336/BioNetGen_2.1.8_rev597.tgz":
        cwd => "/tmp",
        creates => "/tmp/BioNetGen_2.1.8_rev597.tgz",
        require => Package["wget"];
    }
    file { "/tmp/BioNetGen_2.1.8_rev597.tgz": }
    exec { "/bin/tar xzf /tmp/BioNetGen_2.1.8_rev597.tgz":
        cwd => "/home/demo",
        user => demo,
        group => demo,
        creates => "/home/demo/BioNetGen",
        require => [
            File["/tmp/BioNetGen_2.1.8_rev597.tgz"],
            File["/home/demo"],
        ],
    }
    file { "/usr/local/share/BioNetGen":
        ensure => link,
        target => "/home/demo/BioNetGen";
    }
}

class pysb {
#    package {
# TODO
#        "pysb":
#            provider => "pip";
#    }
    user { "demo": ensure => present }
    file { "/home/demo":
        ensure => directory,
        source => "/etc/skel",
        recurse => true,
        owner => "demo",
        group => "demo";
    }
}

class puppet-misc {
    group { "puppet": ensure => present; }
}

class { "apt-update": }
class { "python": }
class { "python-libs": }
class { "perl": }
class { "bionetgen": }
class { "pysb": }
class { "puppet-misc": }

Class["python"] -> Package <| provider == pip |>
Class["apt-update"] -> Package <| provider == apt |>
