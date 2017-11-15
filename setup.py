from distutils.core import setup
setup(
  name = 'pyali',
  packages = ['pyali'], # this must be the same as the name above
  install_requires = ['pandas', 'numpy', 'simplejson'],
  version = '0.1.1',
  description = 'A package for merging alignments',
  author = 'Chris L Tang',
  author_email = 'chris.l.tang@gmail.com',
  url = 'https://github.com/christang/pyali', # use the URL to the github repo
  download_url = 'https://github.com/christang/pyali/archive/0.1.tar.gz', # I'll explain this in a second
  keywords = ['bioinformatics', 'sequence', 'alignment'], # arbitrary keywords
  classifiers = [],
)
