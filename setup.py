import os
from distutils.core import setup
import sys

package_path = os.getcwd()
if not os.path.isdir("./LEBTtunner"):
    os.mkdir("LEBTtunner")
with open("LEBTtunner/package_path.py","w") as f:
    f.write("package_path = '"+package_path+"/LEBTtunner'")


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
