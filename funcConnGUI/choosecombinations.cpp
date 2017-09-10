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

void chooseCombinations::setCheckboxes(int groupNumbers){

    int posY = 1;
    for(int x = 1; x < groupNumbers; ++x){
        for(int y = x+1; y < (groupNumbers+1); ++y){
            QCheckBox *A = new QCheckBox(this);
            A->setText("Group_"+QString::number(x)+"_and_Group_"+QString::number(y));
            A->setChecked(true);
            checkBoxVector.append(A);
            checkBoxVector.last()->setGeometry(20, posY*20, 160, 20);
            posY++;
        }
    }
}


void chooseCombinations::on_buttonBox_accepted()
{
    QVectorIterator<QCheckBox*> i(checkBoxVector);
    //QCheckBox *A = new QCheckBox(this);
    while (i.hasNext()){
        checkedGroups << i.next()->checkState();

    }
}
