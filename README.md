# 課題　 --　RDB10
# planets_data0113

## ①課題内容（どんな作品か）

- 日付と時間を指定すると、該当日時間の天体（水星・金星・太陽・火星・木星・土星・天王星・海王星・冥王星）配置を取得して、その値をjsonファイルとして出力する

https://kotetsu1701.com/blog/msc-planet-position-calculation-1/

上記のサイトに掲載されていたコードを流用（著者の方に連絡して全ての関連ファイルを共有いただきました。感謝！）

  ・lib.py : 天文計算用
  ・m37.py : lib.pyを使用して、任意の天体位置を取得する、その値をjsonファイルで出力させる
  ・from_37pyの中のjsonファイル : m37.pyからの値を取得したjsonファイル。ファイル名には、出力した日付時間が付与される
  
## ②工夫した点・こだわった点

-　参考サイトは、天体位置を時間で値を出力していた※ので、それを角度（360度の円の中では何度なのか）に修正

※天体観測では、時間で出力しても問題ないが、占いチャート（ホロスコープ）上で使用するためには、角度で出力する必要がある

-　参考サイトはではターミナル上で値を確認するだけだったので、値をjsonファイルで出力できるようにした


## ③難しかった点・次回トライしたいこと(又は機能)

- 専門性が高い（エンジニアであり、占いとホロスコープ制作の知見が必要）ので、メンターからの助言が欲しいと思い、下記のサイトでメンターを見つけてアドバイスをお願いしました。
  
https://menta.work/

⚪︎次回は、データベース設計を決める
・暫定的（MVPまで）に考えていること

1. Firebaseを使う
2. 1000人まで登録できる
3. ユーザーは20名の顧客情報まで登録できる
4. 一人の顧客情報は1日1回まで登録できる
5. データは365日間登録され続ける

- [質問と疑問]
   
- [感想]

  chatGPTがコードを出力しても、一部のコードしか出力しないので、そのコードをどこに記載すると正しいループ処理をしてくれるかわからなかった。
  
  そのため、ターミナルで出力された値とjsonで出力された値が違うということが発生。何度かコードの記載箇所を修正して正しい値を出力することができ安心しました。

- [tips]

- [参考記事]

  https://kotetsu1701.com/blog/msc-planet-position-calculation-1/

