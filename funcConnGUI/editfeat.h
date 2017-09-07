#ifndef EDITFEAT_H
#define EDITFEAT_H

#include <QDialog>
#include <QDir>
namespace Ui {
class editFeat;
}

class editFeat : public QDialog
{
    Q_OBJECT

public:
    explicit editFeat(QWidget *parent = 0);
    void setGroups(int groups);
    ~editFeat();

private slots:
    void on_pushButton_Save_clicked();

    void on_comboBox_Groups_currentIndexChanged(int index);

    void on_ld_featDirList_clicked();

private:
    QString WorkingDir = QDir::currentPath();
    QString currentFeatDir;
    QString CommonDir ="/home/";
    Ui::editFeat *ui;
};

#endif // EDITFEAT_H
