#include "choosecombinations.h"
#include "ui_choosecombinations.h"

chooseCombinations::chooseCombinations(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::chooseCombinations)
{
    ui->setupUi(this);
    setWindowTitle(tr("Choose the combinations"));
}

chooseCombinations::~chooseCombinations()
{
    delete ui;
}
