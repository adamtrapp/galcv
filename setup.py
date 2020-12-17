import setuptools

name = 'galcv'

def readme():
    with open('README.md') as f:
        return f.read()

setuptools.setup(name=name,
      version='1.1.1',
      description='A simple calculator for cosmic variance',
      long_description=readme(),
      long_description_content_type="text/markdown",
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
      ],
      keywords='Cosmology Cosmic Variance Galaxies',
      url='https://github.com/adamtrapp/galcv',
      author='Adam Trapp',
      author_email='atrapp@astro.ucla.edu',
      license='MIT',
      packages=setuptools.find_packages(),
      install_requires=['numpy','pandas','scipy'],
      python_requires='~=3.3',
      include_package_data=True,
      zip_safe=False)
