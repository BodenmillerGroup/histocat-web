<template>
  <div>
    <v-snackbar auto-height :color="currentNotificationColor" v-model="show">
      <v-progress-circular class="ma-2" indeterminate v-show="showProgress"></v-progress-circular>
      {{ currentNotificationContent }}
      <v-btn flat @click.native="close">Close</v-btn>
    </v-snackbar>
  </div>
</template>

<script lang="ts">
  import { Component, Vue, Watch } from 'vue-property-decorator';
  import { AppNotification } from '@/modules/main/state';
  import { commitRemoveNotification } from '@/modules/main/mutations';
  import { readFirstNotification } from '@/modules/main/getters';
  import { dispatchRemoveNotification } from '@/modules/main/actions';

  @Component
  export default class NotificationsManager extends Vue {
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
        commitRemoveNotification(this.$store, this.currentNotification);
      }
    }

    get firstNotification() {
      return readFirstNotification(this.$store);
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
          dispatchRemoveNotification(this.$store, { notification: newNotification, timeout: 6500 });
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
