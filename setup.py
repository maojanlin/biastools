from distutils.core import setup
setup(
  name = 'biastools',
  packages = ['biastools'],
  version = '0.0.1',
  license='MIT',
  description = 'The toolkits to analyze reference bias of short DNA read alignment.',
  author = 'Mao-Jan Lin',
  author_email = 'mj.maojanlin@gmail.com',
  url = 'https://github.com/maojanlin/biastools',
  download_url = 'https://github.com/maojanlin/biastools/archive/v_01.tar.gz',
  keywords = ['biastools', 'reference bias', 'alignment'],
  install_requires=[
          'numpy',
          'pysam',
          'pandas',
          'matplotlib',
          'seaborn',
          'scikit-learn',
          'scipy'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
  entry_points={"console_scripts": ["biastools = biastools.biastools:main","biastools_scan = biastools.biastools_scan:main"],},
)

