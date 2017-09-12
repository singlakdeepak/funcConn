#ifndef EDITFEAT_H
#define EDITFEAT_H

#include <QDialog>
#include <QDir>
#include <QList>
namespace Ui {
class editFeat;
}

class editFeat : public QDialog
{
    Q_OBJECT

public:
    explicit editFeat(QWidget *parent = 0);
    void setGroups(int groups);
    QList<QString> get_FileNames();
    ~editFeat();

private slots:
    void on_pushButton_Save_clicked();

    void on_comboBox_Groups_currentIndexChanged(int index);

    void on_ld_featDirList_clicked();

    void on_buttonBox_accepted();


private:
    QString WorkingDir = QDir::currentPath();
    QString currentFeatDir;
    QString CommonDir ="/home/";
    int NGroups;
    QList<QString> FileNames;
    Ui::editFeat *ui;
};

#endif // EDITFEAT_H
