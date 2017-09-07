#include "choosecombinations.h"
#include "ui_choosecombinations.h"
#include<QCheckBox>
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

void chooseCombinations::setCheckboxes(QList<int> groupNumbers){
    QVector<QCheckBox*> checkBoxVector;
    for(int x = 0; x < groupNumbers.at(0); ++x){
        for(int y = 0; y < groupNumbers.at(1); ++y){
            checkBoxVector.append(new QCheckBox(this));
            checkBoxVector.last()->setGeometry(x * 20, y * 20, 20, 20);
        }
    }
}
