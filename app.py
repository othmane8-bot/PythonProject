from flask import Flask, request, render_template, redirect, url_for
from numpy import log as Ln, exp as e

app = Flask(__name__)

# Constantes
CONSTANTS = {
    'V_exp': 1.33e-05,
    'aBA': 194.5302,
    'aAB': -10.7575,
    'lambda_A': 1.127,
    'lambda_B': 0.973,
    'qA': 1.432,
    'qB': 1.4,
    'D_AB': 2.1e-5,
    'D_BA': 2.67e-5
}


def calcul_diffusion(Xa, T):
    if not (0 <= Xa <= 1):
        raise ValueError("La fraction Xa doit être entre 0 et 1")
    if T <= 0:
        raise ValueError("La température doit être positive")

    Xb = 1 - Xa
    phiA = (Xa * CONSTANTS['lambda_A']) / (Xa * CONSTANTS['lambda_A'] + Xb * CONSTANTS['lambda_B'])
    phiB = 1 - phiA
    tauxAB = e(-CONSTANTS['aAB'] / T)
    tauxBA = e(-CONSTANTS['aBA'] / T)
    tetaA = (Xa * CONSTANTS['qA']) / (Xa * CONSTANTS['qA'] + Xb * CONSTANTS['qB'])
    tetaB = 1 - tetaA
    tetaAA = tetaA / (tetaA + tetaB * tauxBA)
    tetaBB = tetaB / (tetaB + tetaA * tauxAB)
    tetaAB = (tetaA * tauxAB) / (tetaA * tauxAB + tetaB)
    tetaBA = (tetaB * tauxBA) / (tetaB * tauxBA + tetaA)

    termes = (
            Xb * Ln(CONSTANTS['D_AB']) +
            Xa * Ln(CONSTANTS['D_BA']) +
            2 * (Xa * Ln(Xa / phiA) + Xb * Ln(Xb / phiB)) +
            2 * Xb * Xa * (
                    (phiA / Xa) * (1 - CONSTANTS['lambda_A'] / CONSTANTS['lambda_B']) +
                    (phiB / Xb) * (1 - CONSTANTS['lambda_B'] / CONSTANTS['lambda_A'])
            ) +
            Xb * CONSTANTS['qA'] * (
                    (1 - tetaBA ** 2) * Ln(tauxBA) +
                    (1 - tetaBB ** 2) * tauxAB * Ln(tauxAB)
            ) +
            Xa * CONSTANTS['qB'] * (
                    (1 - tetaAB ** 2) * Ln(tauxAB) +
                    (1 - tetaAA ** 2) * tauxBA * Ln(tauxBA)
            )
    )
    solution = e(termes)
    erreur = (abs(solution - CONSTANTS['V_exp']) / CONSTANTS['V_exp']) * 100
    return {
        'lnDab': termes,
        'Dab': solution,
        'erreur': round(erreur, 3),
        'Xa': Xa,
        'T': T
    }


@app.route("/")
def home():
    return render_template("home.html")


@app.route("/calcul", methods=["GET", "POST"])
def calcul():
    if request.method == "POST":
        try:
            Xa = float(request.form["Xa"])
            T = float(request.form["T"])
            # Redirect to the result page with query parameters
            return redirect(url_for('result', Xa=Xa, T=T))
        except ValueError as ve:
            return render_template("calcul.html", error=str(ve))
        except Exception as e:
            return render_template("calcul.html", error="Erreur de calcul : " + str(e))
    return render_template("calcul.html")


@app.route("/result")
def result():
    error = None
    try:
        Xa = float(request.args.get("Xa", ""))
        T = float(request.args.get("T", ""))
        data = calcul_diffusion(Xa, T)
    except ValueError as ve:
        error = str(ve)
        data = None
    except Exception as e:
        error = "Erreur de calcul : " + str(e)
        data = None
    return render_template("result.html", result=data, error=error)


@app.route("/explain")
def explain():
    return render_template("explanation.html")


@app.errorhandler(404)
def page_not_found(e):
    return redirect(url_for('home'))


if __name__ == "__main__":
    app.run(debug=True)
