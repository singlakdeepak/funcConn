#ifndef CHOOSECOMBINATIONS_H
#define CHOOSECOMBINATIONS_H

#include <QDialog>

namespace Ui {
class chooseCombinations;
}

class chooseCombinations : public QDialog
{
    Q_OBJECT

public:
    explicit chooseCombinations(QWidget *parent = 0);
    ~chooseCombinations();

private:
    Ui::chooseCombinations *ui;
    //QGroupBox *createNonExclusiveGroup();
};

#endif // CHOOSECOMBINATIONS_H
