<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="headline primary--text">Share Experiment</div>
        <v-spacer/>
        <v-text-field
          v-model="search"
          append-icon="mdi-magnify"
          label="Search"
          single-line
          hide-details
          clearable
        />
      </v-card-title>
      <v-card-text>
        <v-data-table
          :headers="headers"
          :items="users"
          :search="search"
          v-model="selected"
          show-select
        >
        </v-data-table>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn
          @click="submit"
        >
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { Component, Vue } from 'vue-property-decorator';
  import { userModule } from '@/modules/user';
  import { IUserProfile } from '@/modules/user/models';
  import { IShareCreate } from '@/modules/experiment/models';

  @Component
  export default class ShareExperiment extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly userContext = userModule.context(this.$store);

    selected: IUserProfile[] = [];
    search = '';

    headers = [
      {
        text: 'Email',
        sortable: true,
        value: 'email',
        align: 'left',
      },
      {
        text: 'Full Name',
        sortable: true,
        value: 'full_name',
        align: 'left',
      },
    ];

    get users() {
      return this.userContext.getters.adminUsers;
    }

    async mounted() {
      await this.userContext.actions.getUsers();
      this.reset();
    }

    reset() {
      this.selected = [];
    }

    cancel() {
      this.$router.back();
    }

    async submit() {
      const userIds = this.selected.map(item => item.id);
      const data: IShareCreate = {
        user_ids: userIds,
        experiment_id: parseInt(this.$router.currentRoute.params.id, 10)
      };
      await this.experimentContext.actions.createShare(data);
    }
  }
</script>
