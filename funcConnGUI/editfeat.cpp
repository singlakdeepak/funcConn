/* Author: Deepak Singla (singlakdeepak5@gmail.com)
 * Window for selecting FeatDirectories.
 *
 */

#include "editfeat.h"
#include "ui_editfeat.h"
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include<QFileDialog>



editFeat::editFeat(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::editFeat)
{
    ui->setupUi(this);
}

editFeat::~editFeat()
{
    delete ui;
}

void editFeat::setGroups(int groups)
{
    /*
     * It makes the number of groups in the comboBox according to
     * the number of groups added in the main window.
     *
     * INPUT:
     * int groups: Number of groups to be created.
     */
    for (int i =1; i<=groups;i++){
        ui->comboBox_Groups->addItem("Group " + QString::number(i));
    }
}

void editFeat::on_pushButton_Save_clicked()
{
    /*
     * On clicking the save button in the window, it checks
     * whether the folders written row by row in the text editor
     * are present or not. If all files are present, it saves them
     * in a file with the name Group_<IndexNo>_FeatDirs.txt
     *
     */
    int Groupindex = ui->comboBox_Groups->currentIndex();

    QFile file(WorkingDir+"/Group_"+QString::number(Groupindex +1)+"_FeatDirs.txt");
    if (!file.open(QFile::WriteOnly | QFile::Text)){
        // Gives the warning if the following file couldn't be created.
        QMessageBox::warning(this, "title","file not open");
    }
    QTextStream out(&file);
    // Flushes the data extracted from the file into the Plain
    // Text Editor.
    QString text = ui->TextEdit_FeatDirs->toPlainText();
    out << text;
    file.flush();
    file.close();
}

void editFeat::on_comboBox_Groups_currentIndexChanged(int index)
{
    /*
     * When the index on the comboBox is changed, this signal gets activated.
     * It opens the file with name Group_<IndexNo>_FeatDirs.txt in the
     * built application directory if it is present and flushes the contents
     * into the PlainTextEditor.
     *
     */
    ui->TextEdit_FeatDirs->clear();

    QFile file_Feat( WorkingDir+ "/Group_"+QString::number(index+1)+"_FeatDirs.txt");

    // If the following file is present, then only the file is read.
    //Then it flushes the content into Text Editor.
    if (file_Feat.open(QFile::ReadOnly | QFile::Text)){
    QTextStream in_Feat(&file_Feat);
    QString text_Feat = in_Feat.readAll();
    ui->TextEdit_FeatDirs->setPlainText(text_Feat);

    file_Feat.close();
    }
}

void editFeat::on_ld_featDirList_clicked()
{
    /*
     * This signal is for opening a Feat Directories list.
     * The list should be in *.txt format and should contain the
     * absolute paths of the Feat Directories in each row.
     *
     * It checks the String currentFeatDir which stores the
     * initial location from where Feat Directory had been picked.
     * If the current path still exists, then it opens the
     * Directory Chooser from that point else it opens the
     * Directory Chooser from the common path assigned in the scripts.
     */

    QFile Fout(currentFeatDir) ;
    QString new_input;
    if (Fout.exists())
    {

        new_input=QFileDialog::getOpenFileName(this,
                                               tr("Open Feat File"),
                                                    currentFeatDir,
                                                    "Text File (*.txt)"
                                                       );
    }
    else
    {
        new_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open Feat File"),
                                                CommonDir,
                                                "Text File (*.txt)"
                                                   );
    }
    // If 'Cancel' is clicked, then the string will be null,
    // hence it will change the Text Editor data only if some file
    // is chosen.
    if (new_input!=NULL){
        QFile file(new_input);
        if(!file.open(QIODevice::ReadOnly)) {
            QMessageBox::information(0, "error", file.errorString());
        }
        ui->TextEdit_FeatDirs->clear();
        QTextStream in(&file);

        while(!in.atEnd()) {
            QString line = in.readLine();
            ui->TextEdit_FeatDirs->appendPlainText(line);
        }

        file.close();
        currentFeatDir = new_input;
    }
}
