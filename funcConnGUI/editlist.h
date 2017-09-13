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
    QList<QString> get_StructuralFileNames();
    QList<QString> get_FunctionalFileNames();
    ~EditList();


private slots:
    void on_pushButton_Save_clicked();

    void on_comboBox_currentIndexChanged( int index);

    void on_ld_funcFile_clicked();

    void on_ld_StrFile_clicked();

    void on_buttonBox_OK_accepted();

private:
    QString CommonDir ="/home/";
    QString WorkingDir = QDir::currentPath();
    QString CurrentFuncFile,CurrentStructFile;
    QList<QString> StructuralFileNames;
    QList<QString> FunctionalFileNames;
    int NGroups;
    Ui::EditList *ui;
};

#endif // EDITLIST_H
