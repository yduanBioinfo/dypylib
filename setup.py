import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

name = "dypylib"
setuptools.setup(
    name=name, # Replace with your own username
    version="0.0.2",
    author="You Duan",
    author_email="duanyou@outlook.com",
    description="genomics research",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yduanBioinfo/dypylib",
    #packages=setuptools.find_packages(),
    package_dir={name: '.'},
    packages=[name] + ['.'.join((name, x)) for x in setuptools.find_packages()],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
