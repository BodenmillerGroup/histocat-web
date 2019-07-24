<template>
  <div>
    <v-snackbar :color="currentNotificationColor" v-model="show" right>
      <v-progress-circular class="ma-2" indeterminate v-show="showProgress"></v-progress-circular>
      {{ currentNotificationContent }}
      <v-btn text @click.native="close">Close</v-btn>
    </v-snackbar>
  </div>
</template>

<script lang="ts">
  import { mainModule } from '@/modules/main';
  import { AppNotification } from '@/modules/main/models';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  @Component
  export default class NotificationsManager extends Vue {
    readonly mainContext = mainModule.context(this.$store);

    show: boolean = false;
    text: string = '';
    showProgress: boolean = false;
    currentNotification: AppNotification | false = false;

    async hide() {
      this.show = false;
      await new Promise((resolve, reject) => setTimeout(() => resolve(), 500));
    }

    async close() {
      await this.hide();
      await this.removeCurrentNotification();
    }

    async removeCurrentNotification() {
      if (this.currentNotification) {
        this.mainContext.mutations.removeNotification(this.currentNotification);
      }
    }

    get firstNotification() {
      return this.mainContext.getters.firstNotification;
    }

    async setNotification(notification: AppNotification | false) {
      if (this.show) {
        await this.hide();
      }
      if (notification) {
        this.currentNotification = notification;
        this.showProgress = notification.showProgress || false;
        this.show = true;
      } else {
        this.currentNotification = false;
      }
    }

    @Watch('firstNotification')
    async onNotificationChange(
      newNotification: AppNotification | false,
      oldNotification: AppNotification | false,
    ) {
      if (newNotification !== this.currentNotification) {
        await this.setNotification(newNotification);
        if (newNotification) {
          this.mainContext.actions.removeNotification({ notification: newNotification, timeout: 6500 });
        }
      }
    }

    get currentNotificationContent() {
      return this.currentNotification && this.currentNotification.content || '';
    }

    get currentNotificationColor() {
      return this.currentNotification && this.currentNotification.color || 'info';
    }
  }
</script>
