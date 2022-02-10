import os
from distutils.core import setup
import sys

setup(
    name = "LEBTtunner",
    version = "0.0.1",
    author = "Kilean Hwang",
    author_email = "hwang@frib.mse.eud",
    description = ("FRIB LEBT online tunning algorithm"),
    license = "FRIB",
    keywords = "FRIB LEBT",
    url = "",
    packages=['LEBTtunner'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Topic :: Utilities",
        "License :: Free for non-commercial use",
    ],
    zip_safe=False
)
