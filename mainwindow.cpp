#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <math.h>
#include <type_traits>
#include <array>
#include <utility>

template<int N, typename = typename std::enable_if<(N <= 21 && N >= 1)>::type>
struct Factorial
{
    Factorial()
    {
        arr.at(0) = 1L;
        for (unsigned long long i = 1L; i < N; ++i)
        {
            arr.at(i) = arr.at(i - 1L) * i;
        }
    }
    unsigned long long operator()(int n) const
    {
        if (n < 0 || n >= N)
        {
            return 0;
        }
        return arr.at(n);
    }
private:
    std::array<unsigned long long, N> arr;
};

Factorial<21> fact;

template<typename Func>
double Integrate(Func f, double xmin, double xmax, double dx)
{
    double total_area = 0.0;
    double x = xmin;
    auto num_intervals = static_cast<int>(round((xmax - xmin) / dx));
    for(int i = 1; i <= num_intervals; ++i)
    {
        total_area += dx * f(x);
        x += dx;
    }
    return total_area;
}

double LaplasFunction(double x)
{
    if(x == 0.0)
    {
        return 0.0;
    }
    double mul = 1.0;
    if(x < 0)
    {
        x *= -1.0;
        mul = -1.0;
    }
    return 0.3989422804 * Integrate(
        [](double t)
        {
            return exp(-(t*t)/2.0);
        }, 0, x, 0.0001) * mul;
}

inline double GaussFunction(double x)
{
    return 0.3989422804 * exp(-(x*x) / 2.0);
}

int combinations(int n, int m)
{
    if(n < 0 || m < 0)
    {
        return -1;
    }
    if(m > n || n == 0)
    {
        return 0;
    }
    if(n == m || m == 0)
    {
        return 1;
    }

    int n_m = n - m;

    int result = 1;

    if(n_m >= m)
    {
        n_m += 1;
        for(; n_m <= n; ++n_m)
        {
            result *= n_m;
        }
        if(!fact(m))
        {
            return -1;
        }
        return result / fact(m);
    }

    m += 1;
    for(; m <= n; ++m)
    {
        result *= m;
    }
    if(!fact(m))
    {
        return -1;
    }
    return result / fact(n_m);
}

double BernulliFormula(int n, int k, double p)
{
    if(p < 0.0 || p > 1.0)
    {
        return -1.0;
    }
    auto comb = combinations(n, k);
    if(comb < 0)
    {
        return -1.0;
    }
    return comb * pow(p, static_cast<double>(k)) * pow(1.0 - p, static_cast<double>(n - k));
}

double LocalM_L(int n, int k, double p)
{
    if(n < 0 || k < 0 || p < 0.0 || p > 1.0)
    {
        return -1.0;
    }

    double q = 1.0 - p;
    double n_p = n*p;
    double sqrt_npq = sqrt(n_p*q);
    double x = (k - n_p) / sqrt_npq;
    return 0.3989422804 * exp(-(x*x)/2.0) * (1.0 / sqrt_npq);
}

double IntegralM_L(int n, int k1, int k2, double p)
{
    if(n < 0 || k1 < 0 || k2 < 0 || p < 0.0 || p > 1.0)
    {
        return -1.0;
    }

    double q = 1.0 - p;
    double np = n * p;
    double sqrt_npq = sqrt(np*q);

    double x1 = (k1 - np)/sqrt_npq;
    double x2 = (k2 - np)/sqrt_npq;

    return 0.3989422804 * Integrate([](double t)
    {
        return exp(-(t*t) / 2.0);
    }, x1, x2, 0.0001);
}

double PuassonF(int n, int k, double p)
{
    double lambda = n * p;
    if(lambda >= 10.0)
    {
        return -1.0;
    }

    return (pow(lambda, static_cast<double>(k)) * exp(-lambda)) / fact(k);
}

std::pair<bool, double> parseDouble(QString str)
{
    bool success = false;
    str.replace(',', '.');
    double result = str.toDouble(&success);
    return std::make_pair(success, result);
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("TerVerCalculator");

    ui->spinBoxNComb->setMaximum(50);
    ui->spinBoxMComb->setMaximum(50);

    ui->spinBoxNPuasson->setMaximum(1000);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_getAnswerLag_clicked()
{
    QString newText = ") = ";
    ui->labelAnsLag->setText(newText);

    auto parseX = parseDouble(ui->lineEditLag->text());
    if(!parseX.first)
    {
        ui->labelAnsLag->setText(newText.append("<span style=\"color: red;\">error!</span>"));
        return;
    }
    ui->labelAnsLag->setText(newText.append(QString::number(round(LaplasFunction(parseX.second) * 10000)/10000)));
}

void MainWindow::on_getAnswerGauss_clicked()
{
    QString newText = ") = ";
    ui->labelAnsGauss->setText(newText);

    auto parseX = parseDouble(ui->lineEditGauss->text());
    if(!parseX.first)
    {
        ui->labelAnsGauss->setText(newText.append("<span style=\"color: red;\">error!</span>"));
        return;
    }
    ui->labelAnsGauss->setText(newText.append(QString::number(round(GaussFunction(parseX.second) * 10000)/10000)));
}

void MainWindow::on_getAnswerComb_clicked()
{
    QString newText = " = ";
    ui->labelAnsComb->setText(newText);

    int N = ui->spinBoxNComb->value();
    int M = ui->spinBoxMComb->value();
    int ans = combinations(N, M);

    if(ans < 0)
    {
        ui->labelAnsComb->setText(newText.append("<span style=\"color: red;\">error!</span>"));
        return;
    }

    ui->labelAnsComb->setText(newText.append(QString::number(ans)));
}

void MainWindow::on_getAnsBerF_clicked()
{
    QString newText = ") = ";
    ui->labelAnsBerF->setText(newText);
    int n = ui->spinBoxBerFN->value();
    auto parseP = parseDouble(ui->lineEditBerFP->text());
    int k = ui->spinBoxBerFK->value();

    if(!parseP.first || parseP.second < 0.0 || parseP.second > 1.0)
    {
        ui->labelAnsBerF->setText(newText.append("<span style=\"color: red;\">error!</span>"));
        return;
    }

    double ans = BernulliFormula(n, k, parseP.second);
    if(ans < 0.0 || ans > 1.0)
    {
        ui->labelAnsBerF->setText(newText.append("<span style=\"color: red;\">wrong answer!</span>"));
        return;
    }
    newText.append(QString::number(ans));
    ui->labelAnsBerF->setText(newText);
}

void MainWindow::on_getAnsLocalM_L_clicked()
{
    QString newText = ") = ";
    ui->labelAnsLocalM_L->setText(newText);

    int n = ui->spinBoxNLocM_L->value();
    int k = ui->spinBoxKLocM_L->value();
    auto parseP = parseDouble(ui->lineEditPLocM_L->text());

    if(n * parseP.second * (1.0 - parseP.second) < 9.0)
    {
        ui->labelAnsLocalM_L->setText(newText.append("<span style=\"color: red;\">npq &lt; 9</span>"));
        return;
    }

    if(!parseP.first || parseP.second < 0.0 || parseP.second > 1.0)
    {
        ui->labelAnsLocalM_L->setText(newText.append("<span style=\"color: red;\">error!</span>"));
        return;
    }

    double ans = LocalM_L(n, k, parseP.second);

    if(ans < 0.0 || ans > 1.0)
    {
        ui->labelAnsLocalM_L->setText(newText.append("<span style=\"color: red;\">wrong answer!</span>"));
        return;
    }

    newText.append(QString::number(ans));
    ui->labelAnsLocalM_L->setText(newText);
}

void MainWindow::on_getAnsIntM_L_clicked()
{
    QString newText = ") = ";
    ui->labelAnsIntM_L->setText(newText);

    int n = ui->spinBoxNIntM_L->value();
    int k1 = ui->spinBoxK_1IntM_L->value();
    int k2 = ui->spinBoxK_2IntM_L->value();
    auto parseP = parseDouble(ui->lineEditPIntM_L->text());

    if(!parseP.first || parseP.second < 0.0 || parseP.second > 1.0)
    {
        ui->labelAnsIntM_L->setText(newText.append("<span style=\"color: red;\">error!</span>"));
        return;
    }

    if(n * parseP.second * (1.0 - parseP.second) < 9.0)
    {
        ui->labelAnsIntM_L->setText(newText.append("<span style=\"color: red;\">npq &lt; 9</span>"));
        return;
    }

    double ans = IntegralM_L(n, k1, k2, parseP.second);

    if(ans < 0.0 || ans > 1.0)
    {
        ui->labelAnsIntM_L->setText(newText.append("<span style=\"color: red;\">wrong answer!</span>"));
        return;
    }

    newText.append(QString::number(ans));
    ui->labelAnsIntM_L->setText(newText);
}

void MainWindow::on_getAnsPuasson_clicked()
{
    QString newText = ") = ";
    ui->labelAnsPuasson->setText(newText);

    int n = ui->spinBoxNPuasson->value();
    int k = ui->spinBoxKPuasson->value();
    auto parseP = parseDouble(ui->lineEditPPuasson->text());

    if(!parseP.first || parseP.second < 0.0 || parseP.second > 1.0)
    {
        ui->labelAnsPuasson->setText(newText.append("<span style=\"color: red;\">error!</span>"));
        return;
    }

    if(n * parseP.second >= 10.0)
    {
        ui->labelAnsPuasson->setText(newText.append("<span style=\"color: red;\">np &gt;= 10</span>"));
        return;
    }

    double ans = PuassonF(n, k, parseP.second);

    if(ans < 0.0 || ans > 1.0)
    {
        ui->labelAnsPuasson->setText(newText.append("<span style=\"color: red;\">wrong answer!</span>"));
        return;
    }

    newText.append(QString::number(ans));
    ui->labelAnsPuasson->setText(newText);
}
