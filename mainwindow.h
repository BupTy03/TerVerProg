#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_getAnswerLag_clicked();

    void on_getAnswerGauss_clicked();

    void on_getAnswerComb_clicked();

    void on_getAnsBerF_clicked();

    void on_getAnsLocalM_L_clicked();

    void on_getAnsIntM_L_clicked();

    void on_getAnsPuasson_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
