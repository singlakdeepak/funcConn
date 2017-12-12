#ifndef CHOOSECOMBINATIONS_H
#define CHOOSECOMBINATIONS_H

#include <QDialog>
#include <QList>
#include<QCheckBox>
namespace Ui {
class chooseCombinations;
}
//dsada
class chooseCombinations : public QDialog
{
    Q_OBJECT

public:
    explicit chooseCombinations(QWidget *parent = 0);
    ~chooseCombinations();
    void setCheckboxes(int groupNumbers);
    QVector<QCheckBox*> checkBoxVector;
    QList<int> checkedGroups;

private slots:
    void on_buttonBox_accepted();

private:
    Ui::chooseCombinations *ui;
    //QGroupBox *createNonExclusiveGroup();
};

#endif // CHOOSECOMBINATIONS_H
