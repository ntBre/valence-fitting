from flask import Flask
from jinja2 import Environment, PackageLoader, select_autoescape

app = Flask("serve")
env = Environment(
    loader=PackageLoader("serve"), autoescape=select_autoescape()
)


@app.route("/")
def index():
    template = env.get_template("index.html")
    return template.render(
        parameter_ids=["t1", "t2", "t3"],
        molecule_counts=[1, 2, 3],
        pid_counts=[4, 5, 6],
    )
