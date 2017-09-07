#ifndef EDITLIST_H
#define EDITLIST_H

#include <QDialog>
#include <QDir>
namespace Ui {
class EditList;
}

class EditList : public QDialog
{
    Q_OBJECT


public:
    explicit EditList(QWidget *parent = 0);
    void setGroups(int groups);
    ~EditList();


private slots:
    void on_pushButton_Save_clicked();

    void on_comboBox_currentIndexChanged( int index);

    void on_ld_funcFile_clicked();

    void on_ld_StrFile_clicked();

private:
    QString CommonDir ="/home/";
    QString WorkingDir = QDir::currentPath();
    QString CurrentFuncFile,CurrentStructFile;
    Ui::EditList *ui;
};

#endif // EDITLIST_H
