# 使用GitHub Desktop上传代码并协作管理项目的完整操作方案
## 一、前期准备
### 1. 环境配置
- 安装GitHub Desktop：前往[GitHub Desktop官网](https://desktop.github.com/)下载对应系统（Windows/macOS）的客户端并安装。
- 登录GitHub账号：打开GitHub Desktop，通过“File”→“Options”→“Accounts”登录自己的GitHub账号（需提前注册GitHub账号）。
- 确认本地Git配置：GitHub Desktop会自动关联Git，可在“Options”→“Git”中检查用户名和邮箱是否与GitHub账号一致（协作时需保证身份可识别）。

### 2. 远程仓库准备（项目创建者操作）
- 登录GitHub网页端，点击右上角“+”→“New repository”，创建项目仓库：
  - 填写仓库名称（如“team-project”）、描述，选择“Public”/“Private”（私有仓库需邀请协作成员）。
  - 勾选“Add a README file”（初始化仓库，方便协作），点击“Create repository”。
- 邀请协作成员：仓库页面→“Settings”→“Collaborators and teams”→“Add people”，输入成员GitHub账号/邮箱，发送邀请（成员需接受邀请才能协作）。

## 二、核心操作步骤（所有协作成员通用）
### 步骤1：克隆远程仓库到本地
1. 打开GitHub Desktop，点击“File”→“Clone repository”。
2. 在弹出的窗口中：
   - 切换到“GitHub.com”标签，找到创建好的团队仓库（如“team-project”）。
   - 选择本地存储路径（如“D:\Projects\team-project”），点击“Clone”。
   - 等待克隆完成，本地会生成与远程仓库一致的文件夹。

### 步骤2：在本地仓库创建个人专属文件夹
1. 打开克隆后的本地仓库文件夹（如“D:\Projects\team-project”）。
2. 创建以自己用户名/姓名命名的文件夹（如“user-zhangsan”），后续所有个人代码都放在该文件夹内，避免与他人文件冲突。

### 步骤3：本地编写/上传代码到个人文件夹
1. 将自己的代码文件（如.py、.java、.html等）复制到个人专属文件夹中。
2. 回到GitHub Desktop，此时客户端会自动检测到“Changes”（变更）：
   - 左侧“Changes”栏会显示新增的文件夹/文件。
   - 在“Summary”输入框中填写提交说明（如“张三：上传个人模块代码v1.0”），清晰标注作者和内容。
3. 点击“Commit to main”（默认主分支，若需分支开发见进阶操作），将代码提交到本地仓库。

### 步骤4：将本地提交推送到远程仓库
1. GitHub Desktop顶部会显示“Push origin”按钮（若远程仓库有更新，会先提示“Pull origin”，需先拉取再推送）。
2. 点击“Push origin”，等待推送完成，此时个人文件夹和代码会同步到GitHub远程仓库。

### 步骤5：拉取远程仓库最新内容（协作关键）
每次上传代码前，需先拉取远程仓库的最新内容（避免覆盖他人代码）：
1. 打开GitHub Desktop，点击“Pull origin”按钮。
2. 若拉取过程中出现冲突（如多人修改同一文件），客户端会提示“Merge conflicts”，需手动解决冲突后再提交推送。

## 三、进阶操作：分支管理（推荐协作使用）
为避免直接操作主分支导致冲突，建议每人基于主分支创建个人分支：
1. 创建个人分支：
   - GitHub Desktop→“Branch”→“New branch”。
   - 命名分支（如“zhangsan-feature”），选择基于“main”分支创建，点击“Create branch”。
   - 自动切换到个人分支，后续代码提交到该分支。
2. 推送个人分支到远程：完成代码提交后，点击“Push origin”（首次推送分支会自动创建远程分支）。
3. 合并分支到主分支（通过PR）：
   - 登录GitHub网页端→仓库→“Pull requests”→“New pull request”。
   - 选择个人分支（如“zhangsan-feature”）合并到“main”分支，填写PR说明，提交审核。
   - 项目管理员审核通过后，点击“Merge pull request”完成合并，个人代码会同步到主分支。

## 四、协作规范（关键）
1. 每人仅操作自己的专属文件夹，禁止修改他人文件夹内的文件。
2. 每次操作前必须先“Pull origin”拉取最新内容，再“Push origin”推送自己的代码。
3. 提交说明需清晰（格式：【姓名】+ 操作内容 + 版本，如“【李四】修复个人文件夹内登录模块bug v1.1”）。
4. 若需修改公共文件（如README.md、项目配置文件），需先与团队沟通，避免冲突。

## 下一步迭代建议
